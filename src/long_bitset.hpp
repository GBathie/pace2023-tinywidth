#pragma once

#include <bit>
#include <cinttypes>
#include <iterator>
#include <immintrin.h>
#include <array>

template<int N>
class LongBitset
{
private:
	static constexpr uint64_t ONE = 1ULL;
	static constexpr int BIT_SIZE = 64;
	static constexpr int LOG_BITSIZE = 6;
	static constexpr uint64_t LOW_MASK = 0b11'1111;
public:
	static constexpr int MAX_SIZE = 4 * BIT_SIZE * N;
	LongBitset() : b{} {};

	inline bool contains(int i) const { return b[i >> LOG_BITSIZE] & (ONE << (i & LOW_MASK)); }
	inline int  count(int i) const { return contains(i); }
	inline LongBitset &insert(int i) { b[i >> LOG_BITSIZE] |=  (ONE << (i & LOW_MASK)); return *this; }
	inline LongBitset &erase(int i)  { b[i >> LOG_BITSIZE] &= ~(ONE << (i & LOW_MASK)); return *this; }
	inline constexpr int size() const
	{ 
		int res = 0;
		for (auto x: b)
			res += std::popcount(x);
		return  res;
	}
	
	inline bool empty() const
	{
		for (uint64_t x: b)
			if (x)
				return false;
	
		return true;
	}

	static inline LongBitset singleton(int i) { return LongBitset().insert(i); }

	// TODO: OPTIMIZE THIS FUNCTION
	static inline LongBitset full(int n)
	{
		LongBitset res;
		for (int i = 0; i < n; ++i)
			res.insert(i);
		return res;
	}

	LongBitset& operator&=(const LongBitset &other)
	{
		__m256i *ld = (__m256i*)b;
		__m256i const* rd = (__m256i const*) other.b;
		for (int i = 0; i < N; ++i, ++ld, ++rd)
		{
			__m256i l = _mm256_lddqu_si256(ld);
			__m256i r = _mm256_lddqu_si256(rd);
			__m256i res = _mm256_and_si256(l, r);
			_mm256_storeu_si256(ld, res);
		}
		return *this;
	}
	
	LongBitset& operator|=(const LongBitset &other)
	{
		__m256i *ld = (__m256i*)b;
		__m256i const* rd = (__m256i const*) other.b;
		for (int i = 0; i < N; ++i, ++ld, ++rd)
		{
			__m256i l = _mm256_lddqu_si256(ld);
			__m256i r = _mm256_lddqu_si256(rd);
			__m256i res = _mm256_or_si256(l, r);
			_mm256_storeu_si256(ld, res);
		}
		return *this;
	}
	
	LongBitset& operator^=(const LongBitset &other)
	{
		__m256i *ld = (__m256i*)b;
		__m256i const* rd = (__m256i const*) other.b;
		for (int i = 0; i < N; ++i, ++ld, ++rd)
		{
			__m256i l = _mm256_lddqu_si256(ld);
			__m256i r = _mm256_lddqu_si256(rd);
			__m256i res = _mm256_xor_si256(l, r);
			_mm256_storeu_si256(ld, res);	
		}
		return *this;
	}
	
	LongBitset& operator-=(const LongBitset &other)
	{
		__m256i *ld = (__m256i*)b;
		__m256i const* rd = (__m256i const*) other.b;
		for (int i = 0; i < N; ++i, ++ld, ++rd)
		{
			__m256i l = _mm256_lddqu_si256(ld);
			__m256i r = _mm256_lddqu_si256(rd);
			// For some reason, "andnot(x, y)" computes ~x & y
			__m256i res = _mm256_andnot_si256(r, l);
			_mm256_storeu_si256(ld, res);	
		}
		return *this;
	}

	LongBitset operator&(const LongBitset &other) const { LongBitset res = *this; res &= other; return res; }
	LongBitset operator|(const LongBitset &other) const { LongBitset res = *this; res |= other; return res; }
	LongBitset operator^(const LongBitset &other) const { LongBitset res = *this; res ^= other; return res; }
	LongBitset operator-(const LongBitset &other) const { LongBitset res = *this; res -= other; return res; }
	
	LongBitset operator~() const
	{
		LongBitset res = *this;
		for (uint64_t &x: res.b)
			x = ~x;
		return res;
	}

	inline bool operator<=(const LongBitset &other) const
	{
		for (int i = 0; i < 4*N; ++i)
			if ((b[i] & ~other.b[i]) != 0)
				return false;
		return true;
	}

	
	struct Iterator 
    {
    public:
        using iterator_category = std::output_iterator_tag;
        using value_type        = int;
        using reference         = int;

        reference operator*() const { return index; }
        Iterator& operator++() { index = parent.first_set(index + 1); return *this; }  
        Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }
        friend bool operator== (const Iterator& a, const Iterator& b) { return a.index == b.index; };
        friend bool operator!= (const Iterator& a, const Iterator& b) { return a.index != b.index; };  

    private:
    	const LongBitset &parent;
        int index;

    	friend class LongBitset;
    	Iterator(const LongBitset &p, int i) : parent(p), index(i) { }
    };

    Iterator begin() const { return Iterator(*this, first_set()); }
    Iterator begin(int from) const { return Iterator(*this, first_set(from)); }
    Iterator end()   const { return Iterator(*this, 4*N*BIT_SIZE); }
private:
	uint64_t b[4*N];


	// Returns the index of the first bit set to 1,
	// starting from index `from` (inclusive).
	constexpr int first_set(int from = 0) const
	{
		size_t inner = from >> LOG_BITSIZE;
		int outer = from & LOW_MASK;

		while (inner < 4*N)
		{
			int res = first_set_aux(b[inner], outer);
			if (res != BIT_SIZE)
				return res + inner*BIT_SIZE;
			++inner;
			outer = 0;
		}

		return 4*N*BIT_SIZE;
	}

	constexpr int first_set_aux(uint64_t x, int from = 0) const
	{
		uint64_t lb = ONE << from;
		uint64_t mask = lb - 1;

		return std::countr_zero(x & ~mask);
	}
};


#include <iostream>

template<int N>
inline std::ostream &operator<<(std::ostream &os, const LongBitset<N> &b)
{
	os << "b[";
	for (int i: b)
		os << i << ", ";

	return os << "]";
}
