/*
 * Original file written by Farooq Mela (https://github.com/fmela, https://gist.github.com/fmela)
 * Modifications by Sean M. Marks (https://github.com/seanmarks)
 *
 * -----
 *
 * Copyright (c) 2009-2017, Farooq Mela
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// TODO: check for availability of headers
#include <execinfo.h> // for backtrace
#include <dlfcn.h>    // for dladdr
#include <cxxabi.h>   // for __cxa_demangle

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>

namespace debug {

std::string getBacktrace(int skip = 1)
{
	std::ostringstream ss;
	ss << "stack trace:\n";

	void* callstack[128];
	const int max_num_frames = sizeof(callstack) / sizeof(callstack[0]);

	// Get stack addresses
	int num_frames = backtrace(callstack, max_num_frames);
	if (num_frames == 0) {
		ss << "  <empty, possibly corrupt>\n";
		return ss.str();
	}

	// Convert addresses to symbols
	char** symbols = backtrace_symbols(callstack, num_frames);

	char buf[1024];

	try {
		for (int i = skip; i < num_frames; i++) {
			printf("%s\n", symbols[i]);

			Dl_info info;
			if ( dladdr(callstack[i], &info) && info.dli_sname ) {
				// Try to demangle the function name
				char* demangled = NULL;
				int status = -1;
				if (info.dli_sname[0] == '_') {
					demangled = abi::__cxa_demangle(info.dli_sname, NULL, 0, &status);
				}

				snprintf(buf, sizeof(buf), "%-3d %*p %s + %zd\n",
						i,
						int(2 + sizeof(void*) * 2), callstack[i],
						(status == 0 ? demangled :
						               (info.dli_sname == 0 ? symbols[i] : info.dli_sname)),
						(char *)callstack[i] - (char *)info.dli_saddr);

				free(demangled);
			}
			else {
				// Settle for printing the raw output
				snprintf(buf, sizeof(buf), "%-3d %*p %s\n",
						i, int(2 + sizeof(void*) * 2), callstack[i], symbols[i]);
			}

			// Only this line in the try-block can throw
			ss << buf;
		}
	}
	catch (const std::exception& ex) {
		free(symbols);
		std::cerr << "Error in debug::backtrace()";
		throw;
	}

	free(symbols);

	if (num_frames == max_num_frames) {
		ss << "[truncated]\n";
	}

	return ss.str();
}

} // end namespace debug
