// Assert: quick exception throwing with assert-like syntax and informative messages

#ifndef ASSERT_H
#define ASSERT_H

#include <exception>
#include <string>
#include <stdexcept>

#include "utils.h"

#define FANCY_ASSERT(test,message) if (not (test)) { \
  	throw std::runtime_error("assertion failed in " FANCY_FUNCTION "\n" "  " message "\n" "  test: " STRINGIFY(test) "\n"); \
	}

#endif // ifndef ASSERT_H
