/*
 * Copyright (C) 2018, Ibon Guillen (http://giga.cps.unizar.es/~ibon/)
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

#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "bunnykiller.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "Utils/Utils.h"
#include "Filter/Filter.h"
#include "Filter/BoxFilter.h"
#include "Image/RGBColor.h"

/* Forward declarations for ImageIO */
namespace Imaging
{
	template<class T>
	class Image;

	template<class T>
	class ImageRow;

	template<class T>
	void save(const Image<T>& img, const char* filename);
}

namespace Imaging
{
	template<class T = Real>
	class ImageRow
	{
	protected:
		T* m_data;
		size_t m_width;
		size_t m_channels;
	protected:
		friend class Image<T>;

		ImageRow(T* data, size_t w, size_t c) :
			m_data(data),
			m_width(w),
			m_channels(c)
		{
		}
	public:
		inline unsigned width() const
		{
			return m_width;
		}

		inline unsigned channels() const
		{
			return m_channels;
		}
	protected:
		T* data()
		{
			return m_data;
		}
	public:
		inline const T& operator()(unsigned x, unsigned c) const
		{
			return m_data[x * m_channels + c];
		}

		inline T& operator()(unsigned x, unsigned c)
		{
			return m_data[x * m_channels + c];
		}

		const T* data() const
		{
			return m_data;
		}

		size_t size() const
		{
			return m_width * m_channels;
		}
	public:
		/** Add a value to certain pixel in the image row */
		inline void add(unsigned x, const T* values, size_t channels = 1)
		{
			const T* val = &values[0];
			const size_t tchannels = std::min(m_channels, channels);

			std::for_each(&(*this)(x, 0), &(*this)(x, tchannels), [&](T &nth)
			{
				nth += *(val++);
			});
		}

		/** Add a value all pixels in the image row */
		inline void add(const T value)
		{
			std::for_each((*this).data(), (*this).data() + (*this).size(), [&](T& v)
			{
				v += value;
			});
		}

		/** Set a certain pixel in the image row to a value */
		inline void set(unsigned x, const T* values, size_t channels = 1)
		{
			std::copy_n(values, std::min(m_channels, channels), &(*this)(x, 0));
		}

		/** Set all pixels in the image row to a value */
		inline void set(T value)
		{
			std::fill_n((*this).data(), (*this).size(), value);
		}
	public:
		inline void clean()
		{
			set(0.0);
		}
	public:
		void weight(T w);

		void weight(const ImageRow<T>& weigth);
	};

	template<class T>
	void ImageRow<T>::weight(T w)
	{
		std::for_each((*this).data(), (*this).data() + (*this).size(), [&](T& v)
		{
			T tmp = v/w;
			/* Avoid propagating NaNs */
			v = std::isnan(tmp) ? 0.0 : tmp;
		});
	}

	template<class T>
	void ImageRow<T>::weight(const ImageRow<T>& weigth)
	{
		assert((*this).width() == weigth.width());

		if (weigth.channels() == 1) {
			T* vals = (*this).data();

			std::for_each(weigth.data(), weigth.data() + weigth.size(), [&](T w)
			{
				std::for_each(vals, vals + m_channels, [&](T v)
				{
					T tmp = v/w;
					/* Avoid propagating NaNs */
					(*vals++) = std::isnan(tmp) ? 0.0 : tmp;
				});
			});
		} else {
			assert((*this).channels() == weigth.channels());

			T* vals = (*this).data();

			std::for_each(weigth.data(), weigth.data() + weigth.size(), [&](T w)
			{
				T tmp = (*vals)/w;
				/* Avoid propagating NaNs */
				(*vals++) = std::isnan(tmp) ? 0.0 : tmp;
			});
		}
	}

	template<class T = Real>
	class Image
	{
	protected:
		std::unique_ptr<T[]> m_data;
		size_t m_width;
		size_t m_height;
		size_t m_channels;
	public:
		/* Constructors and automatic conversion between images */
		constexpr Image() :
				m_data(),
				m_width(0),
				m_height(0),
				m_channels(0)
		{
		}

		Image(unsigned w, unsigned h, unsigned c = 3) :
				m_width(w),
				m_height(h),
				m_channels(c)
		{
			m_data = std::make_unique<T[]>((*this).size());
		}

		template<class R>
		Image(std::unique_ptr<R[]>&& data, unsigned w, unsigned h, unsigned c) :
				m_width(w),
				m_height(h),
				m_channels(c)
		{
			m_data = std::make_unique<T[]>((*this).size());
			std::copy_n(data.get(), size(), m_data.get());
		}

		Image(std::unique_ptr<T[]>&& data, unsigned w, unsigned h, unsigned c) :
				m_data(std::move(data)),
				m_width(w),
				m_height(h),
				m_channels(c)
		{
		}

		template<class R>
		Image(const Image<R>& img) :
				m_width(img.m_width),
				m_height(img.m_height),
				m_channels(img.m_channels)
		{
			m_data = std::make_unique<T[]>(img.size());
			std::copy_n(img.m_data.get(), img.size(), m_data.get());
		}

		Image(const Image<T>& img) :
				m_width(img.m_width),
				m_height(img.m_height),
				m_channels(img.m_channels)
		{
			m_data = std::make_unique<T[]>(img.size());
			std::copy_n(img.m_data.get(), img.size(), m_data.get());
		}

		template<class R>
		Image(Image<R> && img) :
				m_width(std::move(img.m_width)),
				m_height(std::move(img.m_height)),
				m_channels(std::move(img.m_channels))
		{
			m_data = std::make_unique<T[]>(img.size());
			std::copy_n(img.data(), size(), m_data.get());
		}

		Image(Image<T> && img) :
				m_data(std::move(img.m_data)),
				m_width(std::move(img.m_width)),
				m_height(std::move(img.m_height)),
				m_channels(std::move(img.m_channels))
		{
		}

		template<class R>
		Image& operator=(const Image<R>& img)
		{
			m_data = std::make_unique<T[]>(img.size());
			std::copy_n(img.data(), img.size(), m_data.get());

			m_width = img.m_width;
			m_height = img.m_height;
			m_channels = img.m_channels;

			return *this;
		}

		Image& operator=(const Image<T>& img)
		{
			m_data = std::make_unique<T[]>(img.size());
			std::copy_n(img.data(), img.size(), m_data.get());

			m_width = img.m_width;
			m_height = img.m_height;
			m_channels = img.m_channels;

			return *this;
		}

		template<class R>
		Image& operator=(Image<R> && img)
		{
			m_data = std::make_unique<T[]>(img.size());
			std::copy_n(img.data(), img.size(), m_data.get());

			m_width = img.m_width;
			m_height = img.m_height;
			m_channels = img.m_channels;

			return *this;
		}

		Image& operator=(Image<T> && img)
		{
			m_data = std::move(img.m_data);
			m_width = img.m_width;
			m_height = img.m_height;
			m_channels = img.m_channels;

			return *this;
		}

		~Image()
		{
		}
	public:
		inline unsigned width() const
		{
			return m_width;
		}

		inline unsigned height() const
		{
			return m_height;
		}

		inline unsigned channels() const
		{
			return m_channels;
		}
	protected:
		friend void save<T>(const Image<T>& img, const char* filename);

		T* data()
		{
			return m_data.get();
		}
	public:
		friend class ImageRow<T>;

		ImageRow<T> row(size_t y)
		{
			return ImageRow<T>(&m_data[y * m_width * m_channels], m_width, m_channels);
		}
	public:
		inline const T& operator()(unsigned x, unsigned y, unsigned c) const
		{
			return m_data[y * m_width * m_channels + x * m_channels + c];
		}

		inline T& operator()(unsigned x, unsigned y, unsigned c)
		{
			return m_data[y * m_width * m_channels + x * m_channels + c];
		}

		const T* data() const
		{
			return m_data.get();
		}

		size_t size() const
		{
			return m_width * m_height * m_channels;
		}
	public:
		inline RGBColor operator()(unsigned x, unsigned y) const
		{
			switch (m_channels) {
				case 1:
					return RGBColor(Real((*this)(x, y, 0)));
				case 2:
					return RGBColor(Real((*this)(x, y, 0)), Real((*this)(x, y, 1)), Real(0.0));
				default:
					return RGBColor(Real((*this)(x, y, 0)), Real((*this)(x, y, 1)),
							Real((*this)(x, y, 2)));
			};
		}
	public:
		/** Add a value to certain pixel in the image */
		inline void add(unsigned x, unsigned y, const T* values, size_t channels = 1)
		{
			const T* val = &values[0];
			const size_t tchannels = std::min(m_channels, channels);

			std::for_each(&(*this)(x, y, 0), &(*this)(x, y, tchannels), [&](T &nth)
			{
				nth += *(val++);
			});
		}

		/** Add a value all pixels in the image */
		inline void add(const T value)
		{
			std::for_each((*this).data(), (*this).data() + (*this).size(), [&](T& v)
			{
				v += value;
			});
		}

		/** Set a certain pixel in the image to a value */
		inline void set(unsigned x, unsigned y, const T* values, size_t channels = 1)
		{
			std::copy_n(values, std::min(m_channels, channels), &(*this)(x, y, 0));
		}

		/** Set all pixels in the image to a value */
		inline void set(T value)
		{
			std::fill_n((*this).data(), (*this).size(), value);
		}
	public:
		/** Add a color to certain pixel in the image */
		inline void add(unsigned x, unsigned y, const RGBColor& color)
		{
			add(x, y, &color[0], 3);
		}

		/** Set a certain pixel in the image to a color */
		inline void set(unsigned x, unsigned y, const RGBColor& color)
		{
			set(x, y, &color[0], 3);
		}
	public:
		inline void clean()
		{
			set(0.0);
		}
	public:
		void resize(unsigned w, unsigned h, const Filter* filter = nullptr);

		void transpose(bool reorder_data = true);

		void normalize(T max = 1.);

		void clamp(T min, T max);

		void weight(const Image<T>& weigth);

		void weight(T w);

		Image<T> channels(unsigned c0, unsigned c1) const;
	};

	template<class T>
	void Image<T>::transpose(bool reorder_data)
	{
		if (reorder_data) {
			Image<T> tmp(m_height, m_width, m_channels);

			for (size_t j = 0; j < m_height; j++) {
				for (size_t i = 0; i < m_width; i++) {
					std::copy_n(&(*this)(i, j, 0), m_channels, &tmp(j, i, 0));
				}
			}

			*this = std::move(tmp);
		} else {
			std::swap(m_width, m_height);
		}
	}

	template<class T>
	void Image<T>::normalize(T max)
	{
		size_t size = (*this).size();
		T m_max = *std::max_element(m_data.get(), m_data.get() + size);

		if (m_max > max) {
			std::for_each((*this).data(), (*this).data() + size, [&](T& v)
			{
				v = (v * max) / m_max;
			});
		}
	}

	template<class T>
	void Image<T>::clamp(T min, T max)
	{
		std::for_each((*this).data(), (*this).data() + (*this).size(), [&](T& v)
		{
			v = Utils::clamp(v, min, max);
		});
	}

	template<class T>
	void Image<T>::weight(const Image<T>& weigth)
	{
		assert((*this).width() == weigth.width());
		assert((*this).height() == weigth.height());

		if (weigth.channels() == 1) {
			T* vals = (*this).data();

			std::for_each(weigth.data(), weigth.data() + weigth.size(), [&](T w)
			{
				std::for_each(vals, vals + m_channels, [&](T v)
				{
					T tmp = v/w;
					/* Avoid propagating NaNs */
					(*vals++) = std::isnan(tmp) ? 0.0 : tmp;
				});
			});
		} else {
			assert((*this).channels() == weigth.channels());

			T* vals = (*this).data();

			std::for_each(weigth.data(), weigth.data() + weigth.size(), [&](T w)
			{
				T tmp = (*vals)/w;
				/* Avoid propagating NaNs */
				(*vals++) = std::isnan(tmp) ? 0.0 : tmp;
			});
		}
	}

	template<class T>
	void Image<T>::weight(T w)
	{
		std::for_each((*this).data(), (*this).data() + (*this).size(), [&](T& v)
		{
			T tmp = v/w;
			/* Avoid propagating NaNs */
			v = std::isnan(tmp) ? 0.0 : tmp;
		});
	}

	template<class T>
	Image<T> Image<T>::channels(unsigned c0, unsigned c1) const
	{
		size_t ret_channels = c1 - c0 + 1;
		Image<T> ret(m_width, m_height, ret_channels);

		for (size_t j = 0; j < m_height; j++) {
			for (size_t i = 0; i < m_width; i++) {
				std::copy_n(&(*this)(i, j, c0), ret_channels, &ret(i, j, 0));
			}
		}

		return ret;
	}
}; /* namespace Imaging */

#endif /* _IMAGE_H_ */
