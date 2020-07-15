# AUTHOR: Sean M. Marks (https://github.com/seanmarks)

function(setIfUnset var value)
	if (NOT ${var})
		set(${var} ${value} PARENT_SCOPE)
	endif()
endfunction()
