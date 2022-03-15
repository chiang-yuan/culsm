/* -*- C++ -*------------------------------------------------------------
 *
 *	culsm.h
 *
 * ---------------------------------------------------------------------- */

#ifndef CULSM_H
#define CULSM_H

#define CULSM_VERSION "v2.1 / 15 Mar 2022"

#define WHITESPACE " \t\r\n\f"

#define MAX_NAME	64
#define MAX_LINE	128
#define MAX_STRING	1024
#define MAX_CHUNK	1024

#define BIGINT_FORMAT "%lld"

#define NSECTIONS	25

#define EPS         1e-10

#define BLOCK_SIZE  256

typedef long long bigint;

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#define _USE_MATH_DEFINES

#include <iostream>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>      /* srand, rand */
#include <ctime>        /* time */
#include <cmath>
#include <limits>

#include <cuda.h>
#include <cuda_runtime.h>

#endif // !CULSM_H