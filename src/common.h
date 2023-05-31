#pragma once

#include <iostream>
#include <vector>
#include <list>
#include <unordered_set>
#include <unordered_map>
#include <random>

using contr = std::pair<int,int>;
using contr_seq = std::vector<contr>;
using RNG = std::mt19937;

using RetValue = std::pair<int, contr_seq>;
using MemType = std::unordered_map<std::string, RetValue>;
constexpr int INFTY = 9999;

template <class T, class U>
std::ostream &operator<<(std::ostream &os, const std::pair<T, U> &p)
{
	os << "(" << p.first << ", " << p.second << ")";
	return os;
}


template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v)
{
	os << "[";
	for (const auto &x: v)
		os << x << ", ";
	os << "]";
	return os;
}

template <class T>
std::ostream &operator<<(std::ostream &os, const std::list<T> &v)
{
	os << "{";
	for (const auto &x: v)
		os << x << ", ";
	os << "}";
	return os;
}


template <class T>
std::ostream &operator<<(std::ostream &os, const std::unordered_set<T> &v)
{
	os << "{";
	for (const auto &x: v)
		os << x << ", ";
	os << "}";
	return os;
}

template<class T>
std::unordered_set<T> operator&(const std::unordered_set<T> &l, const std::unordered_set<T> &r)
{
	std::unordered_set<T> res = l;
	for (int x: l)
		if (r.count(x) == 0)
			res.erase(x);

	return res;
}


template<class T>
std::unordered_set<T> operator|(const std::unordered_set<T> &l, const std::unordered_set<T> &r)
{
	std::unordered_set<T> res = l;
	for (int x: r)
		res.insert(x);

	return res;
}


template<class T>
std::unordered_set<T> operator^(const std::unordered_set<T> &l, const std::unordered_set<T> &r)
{
	std::unordered_set<T> res = l | r;
	for (int x: l)
		if (r.count(x) > 0)
			res.erase(x);

	return res;
}

template<class T>
std::unordered_set<T> operator-(const std::unordered_set<T> &l, const std::unordered_set<T> &r)
{
	std::unordered_set<T> res = l;
	for (int x: r)
		res.erase(x);

	return res;
}