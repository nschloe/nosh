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

  std::cout <<
    "Number of vertices: " <<
    mbw->get_number_entities_by_type(0, moab::MBVERTEX) <<
    std::endl;
  std::cout <<
    "Number of edges: " <<
    mbw->get_number_entities_by_type(0, moab::MBEDGE) <<
    std::endl;
  std::cout <<
    "Number of triangles: " <<
    mbw->get_number_entities_by_type(0, moab::MBTRI) <<
    std::endl;
  // get the number of 3D entities
  const int numTets = mbw->get_number_entities_by_type(0, moab::MBTET);
  std::cout << "Number of tetrahedra: " << numTets << std::endl;

  // Create all edges adjacent to tets.
  // Alternative: Create all edges adjacent to vertices.
  if (numTets > 0) {
     const auto tets = mbw->get_entities_by_type(0, moab::MBTET);
     (void) mbw->get_adjacencies(tets, 1, true, moab::Interface::UNION);
  } else {
     const auto tris = mbw->get_entities_by_type(0, moab::MBTRI);
     (void) mbw->get_adjacencies(tris, 1, true, moab::Interface::UNION);
  }
  std::cout <<
    "Number of edges (after creation): " <<
    mbw->get_number_entities_by_type(0, moab::MBEDGE) <<
    std::endl;


#ifndef NDEBUG
  // Print the stats gathered:
  if (0 == global_rank) {
    for (int i = 0; i < nprocs; i++)
      std::cout << " Shared, owned entities on proc " << i << ": " << rbuf[4*i] << " verts, " <<
          rbuf[4*i + 1] << " edges, " << rbuf[4*i + 2] << " faces, " << rbuf[4*i + 3] << " elements" << std::endl;
  }

  // check for illegal vertices
  const auto verts = mbw->get_entities_by_type(0, moab::MBVERTEX);
  for (size_t k = 0; k < verts.size(); k++) {
    const auto edges = mbw->get_adjacencies(
        {verts[k]},
        1,
        false,
        moab::Interface::UNION
        );
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        edges.size() == 0,
        "Node without edge connection. This might be a spurious leftover from " <<
        "the mesh generation and will cause errors when building the system."
        );
  }
#endif

  std::cout << "   mesh_reader >>" << std::endl;
  if (numTets == 0) {
    return std::make_shared<nosh::mesh_tri>(comm, mcomm, mbw->mb);
  }
  return std::make_shared<nosh::mesh_tetra>(comm, mcomm, mbw->mb);
}

}  // namespace nosh
