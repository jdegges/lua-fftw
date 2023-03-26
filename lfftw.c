#define LUA_LIB

#include <lua.h>
#include <lauxlib.h>

#include <fftw3.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LUA_FFTW_NAME "fftw"
#define MT_FFTW_PLAN_1D "fftw plan 1d"


/* Helping to parse/validate arguments passed in on the stack. */

static void
parse_integer (lua_State *L, int index, int required, lua_Integer *var)
{
  int top = lua_gettop (L);

  if (top < index) {
    if (required) {
      lua_pushstring (L, "expected more arguments");
      lua_error (L);
    }
  } else {
    if (!lua_isnumber (L, index)) {
      lua_pushstring (L, "this function requires a number");
      lua_error (L);
    }
    *var = lua_tointeger (L, index);
  }
}

static void
parse_number (lua_State *L, int index, int required, lua_Number *var)
{
  int top = lua_gettop (L);

  if (top < index) {
    if (required) {
      lua_pushstring (L, "expected more arguments");
      lua_error (L);
    }
  } else {
    if (!lua_isnumber (L, index)) {
      lua_pushstring (L, "this function requires a number");
      lua_error (L);
    }
    *var = lua_tonumber (L, index);
  }
}

static void
parse_table (lua_State *L, int index, int required)
{
  int top = lua_gettop (L);

  if (top < index) {
    if (required) {
      lua_pushstring (L, "expected more arguments");
      lua_error (L);
    }
  } else {
    if (!lua_istable (L, index)) {
      lua_pushstring (L, "this function requires a table");
      lua_error (L);
    }
  }
}

static void
parse_userdata (lua_State *L, int index, int required, void **var)
{
  int top = lua_gettop (L);

  if (top < index) {
    if (required) {
      lua_pushstring (L, "expected more arguments");
      lua_error (L);
    }
  } else {
    if (!lua_isuserdata (L, index)) {
      lua_pushstring (L, "this function requires ud");
      lua_error (L);
    }
    *var = lua_touserdata (L, index);
  }
}

/* Wrapper for simple 1D DFT */

struct plan_dft_xd {
  lua_Integer N;
  fftw_plan p;
};

static int
lfftw_plan_dft_1d (lua_State *L)
{
  int n = lua_gettop (L);
  lua_Integer N = 0;
  lua_Integer sign = FFTW_FORWARD;
  lua_Integer flags = FFTW_ESTIMATE;
  fftw_complex* iobuf = NULL;
  fftw_plan p;
  struct plan_dft_xd *ud;

  if (3 < n) {
    lua_pushstring (L, "this function takes at most 3 arguments");
    lua_error (L);
  }

  parse_integer (L, 1, 1, &N);
  parse_integer (L, 2, 0, &sign);
  parse_integer (L, 3, 0, &flags);

  if (0 == N) {
    lua_pushstring (L, "cannot transform an empty buffer");
    lua_error (L);
  }

  iobuf = fftw_malloc (sizeof (fftw_complex) * N);
  if (NULL == iobuf) {
    lua_pushstring (L, "out of memory");
    lua_error (L);
  }

  p = fftw_plan_dft_1d (N, iobuf, iobuf, sign, flags);

  ud = lua_newuserdata (L, sizeof (struct plan_dft_xd));
  if (NULL == ud) {
    lua_pushstring (L, "out of memory");
    lua_error (L);
  }
  ud->N = N;
  ud->p = p;

  /* Associate the metatable with this new userdata object */
  luaL_getmetatable (L, MT_FFTW_PLAN_1D);
  lua_setmetatable (L, -2);

  fftw_free (iobuf);

  return 1;
}

static int
lfftw_execute_dft (lua_State *L)
{
  lua_Integer n = lua_gettop (L);
  lua_Number factor = 1;
  struct plan_dft_xd *ud = NULL;
  lua_Integer N;
  fftw_complex *iobuf;
  fftw_plan p;

  if (4 < n) {
    lua_pushstring (L, "this function takes a maximum of 4 arguments");
    lua_error (L);
  }

  parse_userdata (L, 1, 1, (void **) &ud);
  parse_table (L, 2, 1);
  parse_number (L, 3, 0, &factor);

  if (NULL == ud) {
    lua_pushstring (L, "userdata at argument 1 is NULL?");
    lua_error (L);
  }
  N = ud->N;
  p = ud->p;

  if (N <= 0) {
    lua_pushstring (L, "cannot transform an empty or negative sized buffer");
    lua_error (L);
  }

  if ((size_t) N*2 != lua_rawlen (L, 2)) {
    lua_pushstring (L, "input tables must be the same length");
    lua_error (L);
  }

  iobuf = fftw_malloc (sizeof (fftw_complex) * N);
  if (NULL == iobuf) {
    lua_pushstring (L, "out of memory");
    lua_error (L);
  }

  lua_pushnil (L);
  for (n = 0; 0 != lua_next (L, 2); n++) {
    ((double *)iobuf)[n] = lua_tonumber (L, -1);
    lua_pop (L, 1);
  }

  fftw_execute_dft (p, iobuf, iobuf);

  /* Create new lua array, copy output, and return it to caller */
  lua_createtable (L, N*2, 0);
  for (n = 0; n < N*2; n++) {
    lua_pushinteger (L, n+1);
    lua_pushnumber (L, ((double *)iobuf)[n] * factor);
    lua_settable (L, -3);
  }

  fftw_free (iobuf);

  return 1;
}


static const luaL_Reg fftwlib[] = {
  {"plan_dft_1d", lfftw_plan_dft_1d},
  {NULL, NULL}
};

static const luaL_Reg plan_dft_1d_methods[] = {
  {"execute_dft", lfftw_execute_dft},
  {NULL, NULL}
};

#define SET_FFTW_CONST(s) lua_pushinteger (L, FFTW_##s); lua_setfield (L, -2, #s)

LUALIB_API int luaopen_fftw (lua_State *L)
{
  //luaL_register (L, LUA_FFTW_NAME, fftwlib);
  luaL_newlib(L,fftwlib);
  //lua_newtable(L);
  //luaL_setfuncs(L,fftwlib,0);

  /* Transform sign/direction */
  SET_FFTW_CONST (FORWARD);
  SET_FFTW_CONST (BACKWARD);

  /* Planning-rigor flags */
  SET_FFTW_CONST (ESTIMATE);
  SET_FFTW_CONST (MEASURE);
  SET_FFTW_CONST (PATIENT);
  SET_FFTW_CONST (EXHAUSTIVE);
  SET_FFTW_CONST (WISDOM_ONLY);

  /* Algorithm-restriction flags */
  SET_FFTW_CONST (DESTROY_INPUT);
  SET_FFTW_CONST (PRESERVE_INPUT);
  SET_FFTW_CONST (UNALIGNED);

  if (!luaL_newmetatable (L, MT_FFTW_PLAN_1D)) {
    lua_pushstring (L, "unable to create metatable for simple plans");
    lua_error (L);
  }
  lua_createtable (L, 0, sizeof (plan_dft_1d_methods) / sizeof (luaL_Reg) - 1);
  //luaL_register (L, NULL, plan_dft_1d_methods);
  luaL_setfuncs(L,plan_dft_1d_methods,0);
  lua_setfield (L, -2, "__index");
  lua_pop(L,1); // Pop the metatable so that we can return the library table
  return 1;
}
