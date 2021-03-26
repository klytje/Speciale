function(GLOB_EXCLUDE _output_var _pattern _exclude)
  file(GLOB_RECURSE files ${_pattern})

  foreach (TMP_PATH ${files})
    string (FIND ${TMP_PATH} ${_exclude} EXCLUDE_DIR_FOUND)
    if (NOT ${EXCLUDE_DIR_FOUND} EQUAL -1)
      list (REMOVE_ITEM "files" ${TMP_PATH})
    endif ()
  endforeach(TMP_PATH)
  
  set(${_output_var} ${files} PARENT_SCOPE)
endfunction(GLOB_EXCLUDE)