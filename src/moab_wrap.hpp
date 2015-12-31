// std::vector<double>
// mesh::
// get_data(
//   const std::string & tag_name,
//   const moab::Range & range
//   ) const
// {
//   moab::ErrorCode ierr;
//
//   moab::Tag tag = mbw_->tag_get_handle(tag_name);
//
//   int length;
//   ierr = mbw_->mb->tag_get_length(tag, length);
//   TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);
//
//   moab::DataType type;
//   ierr = mbw_->mb->tag_get_data_type(tag, type);
//   TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);
//
//   TEUCHOS_ASSERT_EQUALITY(type, moab::DataType::MB_TYPE_DOUBLE);
//
//   const int num_data = length * range.size();
//   std::vector<double> data(num_data);
//   ierr = mbw_->mb->tag_get_data(tag, range, &data[0]);
//   TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);
//
//   return data;
// }
#ifndef MOAB_WRAP_HPP
#define MOAB_WRAP_HPP

#include <iterator>
#include <memory>
#include <vector>

#include <moab/Core.hpp>

namespace nosh {
  class moab_wrap {
    public:

      explicit
      moab_wrap(const std::shared_ptr<moab::Core> & _mb):
        mb(_mb)
      {};

      virtual
      ~moab_wrap() {};

      moab::Tag
      tag_get_handle(const std::string & name)
      {
        moab::ErrorCode ierr;
        moab::Tag tag;
        ierr = this->mb->tag_get_handle(name.c_str(), tag);
        TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);
        return tag;
      }

      int
      tag_get_length(const moab::Tag tag)
      {
        moab::ErrorCode ierr;
        int length;
        ierr = this->mb->tag_get_length(tag, length);
        TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);
        return length;
      }

      moab::DataType
      tag_get_data_type(const moab::Tag tag)
      {
        moab::ErrorCode ierr;
        moab::DataType type;
        ierr = this->mb->tag_get_data_type(tag, type);
        TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);
        return type;
      }

      std::vector<double>
      tag_get_data(const moab::Tag tag, const moab::Range & range)
      {
        // TODO assert double
        moab::ErrorCode ierr;
        const int length = this->tag_get_length(tag);
        const int num_data = length * range.size();
        std::vector<double> data(num_data);
        ierr = this->mb->tag_get_data(tag, range, &data[0]);
        TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);
        return data;
      }

      std::vector<int>
      tag_get_int_data(const moab::Tag tag, const moab::Range & range)
      {
        // TODO assert int
        moab::ErrorCode ierr;
        const int length = this->tag_get_length(tag);
        const int num_data = length * range.size();
        std::vector<int> data(num_data);
        ierr = this->mb->tag_get_data(tag, range, &data[0]);
        TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);
        return data;
      }

      std::vector<int>
      tag_get_int_data(
          const moab::Tag tag,
          const std::vector<moab::EntityHandle> & entity_handles
          )
      {
        // TODO assert int
        moab::ErrorCode ierr;
        const int length = this->tag_get_length(tag);
        const int num_data = length * entity_handles.size();
        std::vector<int> data(num_data);
        ierr = this->mb->tag_get_data(
            tag,
            &entity_handles[0],
            entity_handles.size(),
            &data[0]
            );
        TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);
        return data;
      }

      std::vector<double>
      get_coords(const std::vector<moab::EntityHandle> & entity_handles)
      {
        moab::ErrorCode ierr;
        std::vector<double> coords(3 * entity_handles.size());
        ierr = this->mb->get_coords(
            &entity_handles[0],
            entity_handles.size(),
            &coords[0]
            );
        TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS)
        return coords;
      }

      moab::Range
      get_entities_by_dimension(
          const moab::EntityHandle meshset,
          const int dimension,
          const bool recursive=false
          )
      {
        moab::ErrorCode ierr;
        moab::Range entities;
        ierr = this->mb->get_entities_by_dimension(
            meshset,
            dimension,
            entities,
            recursive
            );
        TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);
        return entities;
      }

      moab::Range
      get_adjacencies(
          const moab::Range & from_entities,
          const int to_dimension,
          const bool create_if_missing,
          const int operation_type = moab::Interface::INTERSECT
          )
      {
        moab::ErrorCode ierr;
        moab::Range adj_entities;
        ierr = this->mb->get_adjacencies(
            from_entities,
            to_dimension,
            create_if_missing,
            adj_entities,
            operation_type
            );
        TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);
        return adj_entities;
      }

      moab::Range
      get_adjacencies(
          const std::vector<moab::EntityHandle> & from_entities,
          const int to_dimension,
          const bool create_if_missing,
          const int operation_type = moab::Interface::INTERSECT
          )
      {
        moab::ErrorCode ierr;
        moab::Range adj_entities;
        ierr = this->mb->get_adjacencies(
            &from_entities[0],
            from_entities.size(),
            to_dimension,
            create_if_missing,
            adj_entities,
            operation_type
            );
        TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);
        return adj_entities;
      }

      int
      get_number_entities_by_dimension(
          const moab::EntityHandle meshset,
          const int dimension,
          const bool recursive = false
          )
      {
        moab::ErrorCode ierr;
        int num = 0;
        ierr = this->mb->get_number_entities_by_dimension(
            meshset,
            dimension,
            num,
            recursive
            );
        TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);
        return num;
      }

      std::vector<moab::EntityHandle>
      get_connectivity(
          const moab::EntityHandle entity_handle,
          bool corners_only = false,
          std::vector<moab::EntityHandle> * storage = 0
          )
      {
        moab::ErrorCode ierr;
        const moab::EntityHandle * conn = NULL;
        int numV = 0;
        ierr = this->mb->get_connectivity(
            entity_handle,
            conn,
            numV,
            corners_only,
            storage
            );
        TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);
        // construct vector from array
        std::vector<moab::EntityHandle> conn_vec(numV);
        for (int k = 0; k < numV; k++) {
          conn_vec[k] = conn[k];
        }
        return conn_vec;
      }

    void
    load_file(
        const std::string & file_name,
        const moab::EntityHandle * file_set = 0,
        const std::string & options = "",
        const std::string & set_tag_name = "",
        std::vector<int> set_tag_values = {}
        )
    {
      moab::ErrorCode rval;
      rval = this->mb->load_file(
          file_name.c_str(),
          file_set,
          options.c_str(),
          set_tag_name.c_str(),
          &set_tag_values[0],
          set_tag_values.size()
          );
      TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS);
      return;
    }

    void
    write_mesh(
        const std::string & file_name,
        const moab::EntityHandle * output_list = NULL,
        const int num_sets = 0
        )
    {
      moab::ErrorCode rval;
      rval = this->mb->write_mesh(
          file_name.c_str(),
          output_list,
          num_sets
          );
      TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS);
    }

    public:
      const std::shared_ptr<moab::Core> mb;
  };
}
#endif // MOAB_WRAP_HPP
