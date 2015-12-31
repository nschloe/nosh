// @HEADER
//
//    STK mesh reader.
//    Copyright (C) 2015  Nico Schlömer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// @HEADER
// =============================================================================
#include "mesh_reader.hpp"

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include <moab/Core.hpp>
#include <moab/ParallelComm.hpp>
#include <MBParallelConventions.h>

#include "mesh_tri.hpp"
#include "mesh_tetra.hpp"
#include "moab_wrap.hpp"

namespace nosh
{
std::shared_ptr<nosh::mesh>
read(const std::string & file_name)
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const auto fill_time =
    Teuchos::TimeMonitor::getNewTimer("Nosh: read()");
  Teuchos::TimeMonitor tm(*fill_time);
#endif
  std::cout << ">> mesh_reader" << std::endl;

  const auto comm =
    Teuchos::get_shared_ptr(Teuchos::DefaultComm<int>::getComm());

  const std::string options =
    (comm->getSize() == 1) ?
    "" :
    "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS";

  // Get MOAB instance
  //moab::Interface* mbw->mb = new (std::nothrow) moab::Core;
  auto mb = std::make_shared<moab::Core>();
#ifndef NDEBUG
  TEUCHOS_ASSERT(mb);
#endif
  auto mbw = std::make_shared<moab_wrap>(mb);

  const auto global_rank = comm->getRank();

  //if (nbComms > 1) {
  //  // Split the communicator, into ngroups = nbComms
  //  MPI_Comm_split(MPI_COMM_WORLD, color, global_rank, &comm);
  //}
  //else
  //raw_comm = MPI_COMM_WORLD;
  MPI_Comm raw_comm =
    *(Teuchos::dyn_cast<const Teuchos::MpiComm<int>>(*comm).getRawMpiComm());

  // Get the ParallelComm instance
  auto mcomm = std::make_shared<moab::ParallelComm>(mbw->mb.get(), raw_comm);
  int nprocs = mcomm->proc_config().proc_size();
#ifndef NDEBUG
  MPI_Comm rcomm = mcomm->proc_config().proc_comm();
  assert(rcomm == raw_comm);
#endif

  //MPI_Barrier(MPI_COMM_WORLD);
  MPI_Barrier(raw_comm);

#ifndef NDEBUG
  if (global_rank == 0) {
    int nbComms = 1;
    std::cout << "Reading file " << file_name << "\n with options: " << options <<
         "\n on " << nprocs << " processors on " << nbComms << " communicator(s)\n";
  }
#endif

  // Read the file with the specified options
  mbw->load_file(file_name, 0, options);

  moab::ErrorCode rval;
  moab::Range shared_ents;
  // Get entities shared with all other processors
  rval = mcomm->get_shared_entities(-1, shared_ents);
  TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS)

  // Filter shared entities with not not_owned, which means owned
  moab::Range owned_entities;
  rval = mcomm->filter_pstatus(
      shared_ents,
      PSTATUS_NOT_OWNED,
      PSTATUS_NOT,
      -1,
      &owned_entities
      );
  TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS)

  unsigned int nums[4] = {0}; // to store the owned entities per dimension
  for (int i = 0; i < 4; i++)
    nums[i] = (int)owned_entities.num_of_dimension(i);
  std::vector<int> rbuf(nprocs*4, 0);
  MPI_Gather(nums, 4, MPI_INT, &rbuf[0], 4, MPI_INT, 0, raw_comm);

#ifndef NDEBUG
  // Print the stats gathered:
  if (0 == global_rank) {
    for (int i = 0; i < nprocs; i++)
      std::cout << " Shared, owned entities on proc " << i << ": " << rbuf[4*i] << " verts, " <<
          rbuf[4*i + 1] << " edges, " << rbuf[4*i + 2] << " faces, " << rbuf[4*i + 3] << " elements" << std::endl;
  }
#endif

  // get the number of 3D entities
  const int num3d = mbw->get_number_entities_by_dimension(0, 3);

  std::cout << "   mesh_reader >>" << std::endl;
  if (num3d == 0) {
    return std::make_shared<nosh::mesh_tri>(comm, mcomm, mbw->mb);
  }
  return std::make_shared<nosh::mesh_tetra>(comm, mcomm, mbw->mb);
}

}  // namespace nosh
