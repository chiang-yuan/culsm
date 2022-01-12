#include "fix.h"

Fix::Fix(Error *error_)
{
	error = error_;
}

Fix::Fix(Error *error_, int type_, double dispx_, double dispy_, double dispz_)
{
	error = error_;

	type = type_;
	dispx = dispx_;
	dispy = dispy_;
	dispz = dispz_;
}

Fix::~Fix()
{
}