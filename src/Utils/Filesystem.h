/*
 * Copyright (C) 2017, Ibon Guillen (http://giga.cps.unizar.es/~ibon/)
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom
 * the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _UTILS_FILESYSTEM_H_
#define _UTILS_FILESYSTEM_H_

#include <cstring>
#include <string>
#if defined(WIN32)
#include <direct.h>
#include <regex>
#else
#include <sys/types.h>
#include <sys/stat.h>
#endif

namespace Filesystem
{
	/* Returns the extension of a given filename */
	std::string file_extension(const char* filename)
	{
		const char* extension = strrchr(filename, '.');
		if (extension != nullptr) {
			return std::string(extension+1);
		}

		return std::string("");
	}

	/* Returns the root path of a given filename */
	std::string file_path(const char* filename)
	{
		const char* path = strrchr(filename, '/');
#if defined(WIN32)
		if (path == nullptr) {
			path = strrchr(filename, '\\');
		}
#endif
		if (path != nullptr) {
			return std::string(filename, size_t(path - filename));
		}

		return std::string("");
	}

	/* Returns the base name of a given filename */
	std::string base_name(const char* filename, bool include_path = false)
	{
		const char* path = nullptr;
		if (!include_path) {
			path = strrchr(filename, '/');
#if defined(WIN32)
			if (path == nullptr) {
				path = strrchr(filename, '\\');
			}
#endif
			/* Omit the file separator */
			path = path+1;
		}

		/* If there is no path, use the entire filename */
		path = (path == nullptr) ? filename : path;

		const char* extension = strrchr(filename, '.');
		if (extension == nullptr) {
			return std::string(path);
		} else {
			return std::string(path, size_t(extension - path));
		}
	}

	/* Creates the folder structure of a given path (if it doesn't exist) */
	void create_directory(const char* filename)
	{
		/* Create directory structure recursively, starting from root */
		std::string root = Filesystem::file_path(filename);
		if (!root.empty()) {
			create_directory(root.c_str());
		}
#if defined(WIN32)
		/* Avoid device roots */
		if (std::regex_match(filename, std::regex("[A-Z]:"))) {
			return;
		}

		int res = _mkdir(filename);
#else
		mode_t mode = 0775;
		int res = mkdir(filename, mode);
#endif
		/* Check errors and ignore if directory exists */
		int err = errno;

		if (res != 0 && err != EEXIST) {
			throw std::runtime_error("Error creating directory " + std::string(filename));
		}
	}
}

#endif /* _UTILS_FILESYSTEM_H_ */
