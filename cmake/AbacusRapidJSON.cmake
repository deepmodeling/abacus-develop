function(abacus_configure_rapidjson target_name)
  find_package(RapidJSON CONFIG REQUIRED)

  if(TARGET RapidJSON::RapidJSON)
    target_link_libraries("${target_name}" INTERFACE RapidJSON::RapidJSON)
  elseif(TARGET rapidjson)
    target_link_libraries("${target_name}" INTERFACE rapidjson)
  elseif(TARGET RapidJSON)
    target_link_libraries("${target_name}" INTERFACE RapidJSON)
  elseif(DEFINED RapidJSON_INCLUDE_DIRS AND NOT "${RapidJSON_INCLUDE_DIRS}" STREQUAL "")
    target_include_directories("${target_name}" INTERFACE ${RapidJSON_INCLUDE_DIRS})
  elseif(DEFINED RapidJSON_INCLUDE_DIR AND NOT "${RapidJSON_INCLUDE_DIR}" STREQUAL "")
    target_include_directories("${target_name}" INTERFACE "${RapidJSON_INCLUDE_DIR}")
  else()
    message(
      FATAL_ERROR
      "RapidJSON was found, but it did not provide a usable target or include directory. "
      "Expected one of RapidJSON::RapidJSON, rapidjson, RapidJSON, RapidJSON_INCLUDE_DIRS, or RapidJSON_INCLUDE_DIR."
    )
  endif()

  abacus_add_feature_definitions(__RAPIDJSON)
endfunction()
