// std::vector<double>
// mesh::
// get_data(
//   const std::string & tag_name,
//   const moab::Range & range
//   ) const
// {
//   moab::ErrorCode rval;
//
//   moab::Tag tag = mbw_->tag_get_handle(tag_name);
//
//   int length;
//   rval = mbw_->mb->tag_get_length(tag, length);
//   TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS);
//
//   moab::DataType type;
//   rval = mbw_->mb->tag_get_data_type(tag, type);
//   TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS);
//
//   TEUCHOS_ASSERT_EQUALITY(type, moab::DataType::MB_TYPE_DOUBLE);
//
//   const int num_data = length * range.size();
//   std::vector<double> data(num_data);
//   rval = mbw_->mb->tag_get_data(tag, range, &data[0]);
//   TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS);
//
//   return data;
// }
#ifndef MOAB_WRAP_HPP
#define MOAB_WRAP_HPP

#include <iterator>
#include <memory>
#include <sstream>
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
        moab::Tag tag;
        const auto rval = this->mb->tag_get_handle(name.c_str(), tag);
        if (rval != moab::MB_SUCCESS) {
          std::ostringstream oss;
          oss << "error in moab::tag_get_handle(string) "
              << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
          throw std::runtime_error(oss.str());
        }
        return tag;
      }

      int
      tag_get_length(const moab::Tag tag)
      {
        int length;
        const auto rval = this->mb->tag_get_length(tag, length);
        if (rval != moab::MB_SUCCESS) {
          std::ostringstream oss;
          oss << "error in moab::tag_get_length "
              << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
          throw std::runtime_error(oss.str());
        }
        return length;
      }

      moab::DataType
      tag_get_data_type(const moab::Tag tag)
      {
        moab::DataType type;
        const auto rval = this->mb->tag_get_data_type(tag, type);
        if (rval != moab::MB_SUCCESS) {
          std::ostringstream oss;
          oss << "error in moab::tag_get_data_type "
              << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
          throw std::runtime_error(oss.str());
        }
        return type;
      }

      std::vector<double>
      tag_get_data(
          const moab::Tag tag,
          const moab::Range & range
          )
      {
        // TODO assert double
        const int length = this->tag_get_length(tag);
        const int num_data = length * range.size();
        std::vector<double> data(num_data);
        const auto rval = this->mb->tag_get_data(tag, range, data.data());
        if (rval != moab::MB_SUCCESS) {
          std::ostringstream oss;
          oss << "error in moab::tag_get_data "
              << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
          throw std::runtime_error(oss.str());
        }
        return data;
      }

      std::vector<int>
      tag_get_int_data(
          const moab::Tag tag,
          const moab::Range & range
          )
      {
        // TODO assert int
        const int length = this->tag_get_length(tag);
        const int num_data = length * range.size();
        std::vector<int> data(num_data);
        const auto rval = this->mb->tag_get_data(tag, range, data.data());
        if (rval != moab::MB_SUCCESS) {
          std::ostringstream oss;
          oss << "error in moab::tag_get_int_data "
              << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
          throw std::runtime_error(oss.str());
        }
        return data;
      }

      std::vector<int>
      tag_get_int_data(
          const moab::Tag tag,
          const std::vector<moab::EntityHandle> & entity_handles
          )
      {
        // TODO assert int
        const int length = this->tag_get_length(tag);
        const int num_data = length * entity_handles.size();
        std::vector<int> data(num_data);
        const auto rval = this->mb->tag_get_data(
            tag,
            entity_handles.data(),
            entity_handles.size(),
            data.data()
            );
        if (rval != moab::MB_SUCCESS) {
          std::ostringstream oss;
          oss << "error in moab::tag_get_int_data "
              << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
          throw std::runtime_error(oss.str());
        }
        return data;
      }

      std::vector<double>
      get_coords(const std::vector<moab::EntityHandle> & entity_handles)
      {
        std::vector<double> coords(3 * entity_handles.size());
        const auto rval = this->mb->get_coords(
            entity_handles.data(),
            entity_handles.size(),
            coords.data()
            );
        if (rval != moab::MB_SUCCESS) {
          std::ostringstream oss;
          oss << "error in moab::get_coords "
              << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
          throw std::runtime_error(oss.str());
        }
        return coords;
      }

      moab::Range
      get_entities_by_type(
          const moab::EntityHandle meshset,
          const moab::EntityType type,
          const bool recursive=false
          )
      {
        moab::Range entities;
        const auto rval = this->mb->get_entities_by_type(
            meshset,
            type,
            entities,
            recursive
            );
        if (rval != moab::MB_SUCCESS) {
          std::ostringstream oss;
          oss << "error in moab::get_entities_by_type "
              << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
          throw std::runtime_error(oss.str());
        }
        return entities;
      }

      moab::Range
      get_entities_by_dimension(
          const moab::EntityHandle meshset,
          const int dimension,
          const bool recursive=false
          )
      {
        moab::Range entities;
        const auto rval = this->mb->get_entities_by_dimension(
            meshset,
            dimension,
            entities,
            recursive
            );
        if (rval != moab::MB_SUCCESS) {
          std::ostringstream oss;
          oss << "error in moab::get_entities_by_dimension "
              << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
          throw std::runtime_error(oss.str());
        }
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
        moab::Range adj_entities;
        const auto rval = this->mb->get_adjacencies(
            from_entities,
            to_dimension,
            create_if_missing,
            adj_entities,
            operation_type
            );
        if (rval != moab::MB_SUCCESS) {
          std::ostringstream oss;
          oss << "error in moab::get_adjacencies "
              << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
          throw std::runtime_error(oss.str());
        }
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
        moab::Range adj_entities;
        const auto rval = this->mb->get_adjacencies(
            from_entities.data(),
            from_entities.size(),
            to_dimension,
            create_if_missing,
            adj_entities,
            operation_type
            );
        if (rval != moab::MB_SUCCESS) {
          std::ostringstream oss;
          oss << "error in moab::get_adjacencies "
              << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
          throw std::runtime_error(oss.str());
        }
        return adj_entities;
      }

      int
      get_number_entities_by_dimension(
          const moab::EntityHandle meshset,
          const int dimension,
          const bool recursive = false
          )
      {
        int num = 0;
        const auto rval = this->mb->get_number_entities_by_dimension(
            meshset,
            dimension,
            num,
            recursive
            );
        if (rval != moab::MB_SUCCESS) {
          std::ostringstream oss;
          oss << "error in moab::get_number_entities_by_dimension "
              << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
          throw std::runtime_error(oss.str());
        }
        return num;
      }

      int
      get_number_entities_by_type(
          const moab::EntityHandle meshset,
          const moab::EntityType type,
          const bool recursive = false
          )
      {
        int num = 0;
        const auto rval = this->mb->get_number_entities_by_type(
            meshset,
            type,
            num,
            recursive
            );
        if (rval != moab::MB_SUCCESS) {
          std::ostringstream oss;
          oss << "error in moab::get_number_entities_by_type "
              << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
          throw std::runtime_error(oss.str());
        }
        return num;
      }

      std::vector<moab::EntityHandle>
      get_connectivity(
          const moab::EntityHandle entity_handle,
          bool corners_only = false,
          std::vector<moab::EntityHandle> * storage = 0
          )
      {
        const moab::EntityHandle * conn = NULL;
        int numV = 0;
        const auto rval = this->mb->get_connectivity(
            entity_handle,
            conn,
            numV,
            corners_only,
            storage
            );
        if (rval != moab::MB_SUCCESS) {
          std::ostringstream oss;
          oss << "error in moab::get_connectivity "
              << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
          throw std::runtime_error(oss.str());
        }
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
        const moab::EntityHandle * file_set = nullptr,
        const std::string & options = "",
        const std::string & set_tag_name = "",
        std::vector<int> set_tag_values = {}
        )
    {
      const auto rval = this->mb->load_file(
          file_name.c_str(),
          file_set,
          options.c_str(),
          set_tag_name.c_str(),
          set_tag_values.data(),
          set_tag_values.size()
          );
      if (rval != moab::MB_SUCCESS) {
        std::ostringstream oss;
        oss << "error in moab::load_file "
            << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
        throw std::runtime_error(oss.str());
      }
      return;
    }

    void
    write_file(
        const std::string & file_name,
        const std::string & file_type = "",
        const std::string & options = "",
        const std::vector<moab::EntityHandle> & output_sets = {},
        const std::vector<moab::Tag> & tag_list = {}
        )
    {
      const auto rval = this->mb->write_file(
          file_name.c_str(),
          file_type.c_str(),
          options.c_str(),
          output_sets.data(),
          output_sets.size(),
          tag_list.data(),
          tag_list.size()
          );
      if (rval != moab::MB_SUCCESS) {
        std::ostringstream oss;
        oss << "error in moab::write_file "
            << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
        throw std::runtime_error(oss.str());
      }
    }

    std::tuple<moab::Tag, bool>
    tag_get_handle(
        const std::string & name,
        int size,
        moab::DataType type,
        unsigned flags = 0,
        const void * default_value = 0
        )
    {
      moab::Tag tag_handle;
      bool created;
      const auto rval = this->mb->tag_get_handle(
          name.c_str(),
          size,
          type,
          tag_handle,
          flags,
          default_value,
          &created
          );
      if (rval != moab::MB_SUCCESS) {
        std::ostringstream oss;
        oss << "error in moab::tag_get_handle(string, size, ...) "
            << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
        throw std::runtime_error(oss.str());
      }
      return std::make_tuple(tag_handle, created);
    }

    void
    tag_set_data(
        moab::Tag tag_handle,
        const std::vector<moab::EntityHandle> & entity_handles,
        const void * tag_data
        )
    {
      const auto rval = this->mb->tag_set_data(
        tag_handle,
        &entity_handles[0],
        entity_handles.size(),
        tag_data
        );
      if (rval != moab::MB_SUCCESS) {
        std::ostringstream oss;
        oss << "error in moab::tag_set_data "
            << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
        throw std::runtime_error(oss.str());
      }
      return;
    }

    void
    tag_set_data(
        moab::Tag tag_handle,
        const moab::Range & entity_handles,
        const void * tag_data
        )
    {
      const auto rval = this->mb->tag_set_data(
        tag_handle,
        entity_handles,
        tag_data
        );
      if (rval != moab::MB_SUCCESS) {
        std::ostringstream oss;
        oss << "error in moab::tag_set_data "
            << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
        throw std::runtime_error(oss.str());
      }
      return;
    }

    moab::EntityHandle
    create_meshset(
        const unsigned int options,
        int start_id = 0
        )
    {
      moab::EntityHandle handle;
      const auto rval = this->mb->create_meshset(
          options,
          handle,
          start_id
          );
      if (rval != moab::MB_SUCCESS) {
        std::ostringstream oss;
        oss << "error in moab::create_meshset "
            << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
        throw std::runtime_error(oss.str());
      }
      return handle;
    }

    void
    add_entities(
        const moab::EntityHandle & meshset,
        const moab::Range & entities
        )
    {
      const auto rval = this->mb->add_entities(
          meshset,
          entities
          );
      if (rval != moab::MB_SUCCESS) {
        std::ostringstream oss;
        oss << "error in moab::add_entities "
            << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
        throw std::runtime_error(oss.str());
      }
      return;
    }

    void
    add_entities(
        const moab::EntityHandle & meshset,
        const std::vector<moab::EntityHandle> & entities
        )
    {
      const auto rval = this->mb->add_entities(
          meshset,
          &entities[0],
          entities.size()
          );
      if (rval != moab::MB_SUCCESS) {
        std::ostringstream oss;
        oss << "error in moab::add_entities "
            << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
        throw std::runtime_error(oss.str());
      }
      return;
    }

    void
    unite_meshset(
        moab::EntityHandle meshset1,
        const moab::EntityHandle meshset2
        )
    {
      const auto rval = this->mb->unite_meshset(
          meshset1,
          meshset2
          );
      if (rval != moab::MB_SUCCESS) {
        std::ostringstream oss;
        oss << "error in moab::unite_meshset "
            << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
        throw std::runtime_error(oss.str());
      }
      return;
    }

    int
    get_dimension()
    {
      int dim;
      const auto rval = this->mb->get_dimension(dim);
      if (rval != moab::MB_SUCCESS) {
        std::ostringstream oss;
        oss << "error in moab::get_dimension "
            << "(error code " << rval << ", " << translate_error_code_(rval) << ")";
        throw std::runtime_error(oss.str());
      }
      return dim;
    }

    bool
    contains_entities(
        const moab::EntityHandle & meshset,
        const std::vector<moab::EntityHandle> & entities,
        const int operation_type = moab::Interface::INTERSECT
        )
    {
      return this->mb->contains_entities(
          meshset,
          &entities[0],
          entities.size(),
          operation_type
          );
    }

    public:
      const std::shared_ptr<moab::Core> mb;

    private:
      std::string
      translate_error_code_(const moab::ErrorCode error_code)
      {
        switch (error_code) {
          case moab::MB_INDEX_OUT_OF_RANGE:
            return "MB_INDEX_OUT_OF_RANGE";
          case moab::MB_TYPE_OUT_OF_RANGE:
            return "MB_TYPE_OUT_OF_RANGE";
          case moab::MB_MEMORY_ALLOCATION_FAILED:
            return "MB_MEMORY_ALLOCATION_FAILED";
          case moab::MB_ENTITY_NOT_FOUND:
            return "MB_ENTITY_NOT_FOUND";
          case moab::MB_MULTIPLE_ENTITIES_FOUND:
            return "MB_MULTIPLE_ENTITIES_FOUND";
          case moab::MB_TAG_NOT_FOUND:
            return "MB_TAG_NOT_FOUND";
          case moab::MB_FILE_DOES_NOT_EXIST:
            return "MB_FILE_DOES_NOT_EXIST";
          case moab::MB_FILE_WRITE_ERROR:
            return "MB_FILE_WRITE_ERROR";
          case moab::MB_NOT_IMPLEMENTED:
            return "MB_NOT_IMPLEMENTED";
          case moab::MB_ALREADY_ALLOCATED:
            return "MB_ALREADY_ALLOCATED";
          case moab::MB_VARIABLE_DATA_LENGTH:
            return "MB_VARIABLE_DATA_LENGTH";
          case moab::MB_INVALID_SIZE:
            return "MB_INVALID_SIZE";
          case moab::MB_UNSUPPORTED_OPERATION:
            return "MB_UNSUPPORTED_OPERATION";
          case moab::MB_UNHANDLED_OPTION:
            return "MB_UNHANDLED_OPTION";
          case moab::MB_STRUCTURED_MESH:
            return "MB_STRUCTURED_MESH";
          case moab::MB_FAILURE:
            return "MB_FAILURE";
          default:
            return "unknown error";
        }
      }
  };
}
#endif // MOAB_WRAP_HPP
