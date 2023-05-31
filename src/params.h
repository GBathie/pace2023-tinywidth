#pragma once

// Params
constexpr int LB_K = 25;
#ifdef OPTIL
	#pragma message("Compiling OPTIL version")
	constexpr int LB_TIME_S = 60;
#else
	constexpr int LB_TIME_S = 5;
#endif

constexpr int DEFAULT_OUTER_SP = 10;
constexpr int DEFAULT_INNER_SP = 50;

constexpr int DEFAULT_TREE = 200;
constexpr int DEFAULT_OUTER = 100;
constexpr int DEFAULT_INNER = 200;
