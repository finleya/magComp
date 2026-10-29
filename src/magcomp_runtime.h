#ifndef MAGCOMP_RUNTIME_H
#define MAGCOMP_RUNTIME_H

#include <cstdio>
#include <exception>
#include <stdexcept>
#include <type_traits>
#include <Rinternals.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>

namespace magcomp {

// R errors/interrupts normally longjmp past C++ destructors. Protect calls
// into R below the C++ frames and resume the R condition only at the boundary.
struct RUnwind {};
inline SEXP continuation = NULL;

inline void unwind(void *, Rboolean jump)
{
  if(jump) throw RUnwind();
}

template<class F> auto rCall(F function)
{
  using Result = decltype(function());
  if constexpr(std::is_same_v<Result, SEXP>){
    return R_UnwindProtect([](void *data) -> SEXP {
      return (*static_cast<F *>(data))();
    }, &function, unwind, NULL, continuation);
  }else if constexpr(std::is_void_v<Result>){
    R_UnwindProtect([](void *data) -> SEXP {
      (*static_cast<F *>(data))();
      return R_NilValue;
    }, &function, unwind, NULL, continuation);
  }else{
    struct Call { F &function; Result result; } call{function, {}};
    R_UnwindProtect([](void *data) -> SEXP {
      Call &call = *static_cast<Call *>(data);
      call.result = call.function();
      return R_NilValue;
    }, &call, unwind, NULL, continuation);
    return call.result;
  }
}

template<class... Args>
[[noreturn]] void error(const char *format, Args... args)
{
  char message[2048];
  if constexpr(sizeof...(args) == 0)
    std::snprintf(message, sizeof(message), "%s", format);
  else
    std::snprintf(message, sizeof(message), format, args...);
  throw std::runtime_error(message);
}

template<class F> SEXP nativeCall(F function)
{
  SEXP token = PROTECT(R_MakeUnwindCont());
  SEXP previous = continuation;
  continuation = token;
  SEXP result = R_NilValue;
  bool interrupted = false;
  char message[2048] = {};
  try {
    result = function();
  }catch(const RUnwind &) {
    interrupted = true;
  }catch(const std::exception &exception) {
    std::snprintf(message, sizeof(message), "%s", exception.what());
  }catch(...) {
    std::snprintf(message, sizeof(message), "unknown native exception");
  }
  continuation = previous;
  if(interrupted) R_ContinueUnwind(token);
  UNPROTECT(1);
  if(message[0]) Rf_error("%s", message);
  return result;
}

// Keep RNG cleanup outside the C++ frames: PutRNGstate can itself allocate.
// Implementations set rngActive only after GetRNGstate succeeds.
template<class F> SEXP nativeCallWithRng(F function)
{
  struct Call { F &function; bool rngActive; } call{function, false};
  return R_UnwindProtect([](void *data) -> SEXP {
    Call &call = *static_cast<Call *>(data);
    return nativeCall([&]() { return call.function(call.rngActive); });
  }, &call, [](void *data, Rboolean) {
    Call &call = *static_cast<Call *>(data);
    if(call.rngActive){
      call.rngActive = false;
      PutRNGstate();
    }
  }, &call, NULL);
}

// These wrappers also cover allocation failures while native vectors exist,
// including materialization of ALTREP inputs.
#define MAGCOMP_R_WRAPPER(name) \
  template<class... Args> auto guarded_##name(Args... args) { \
    return rCall([&]() { return name(args...); }); \
  }
MAGCOMP_R_WRAPPER(Rf_allocVector)
MAGCOMP_R_WRAPPER(Rf_allocMatrix)
MAGCOMP_R_WRAPPER(Rf_allocArray)
MAGCOMP_R_WRAPPER(Rf_coerceVector)
MAGCOMP_R_WRAPPER(Rf_setAttrib)
MAGCOMP_R_WRAPPER(Rf_mkChar)
MAGCOMP_R_WRAPPER(Rf_mkString)
MAGCOMP_R_WRAPPER(Rf_ScalarReal)
MAGCOMP_R_WRAPPER(Rf_ScalarInteger)
MAGCOMP_R_WRAPPER(Rf_asReal)
MAGCOMP_R_WRAPPER(Rf_asInteger)
MAGCOMP_R_WRAPPER(Rf_asLogical)
MAGCOMP_R_WRAPPER(Rf_asChar)
MAGCOMP_R_WRAPPER(Rf_getAttrib)
MAGCOMP_R_WRAPPER(Rf_length)
MAGCOMP_R_WRAPPER(REAL)
MAGCOMP_R_WRAPPER(INTEGER)
MAGCOMP_R_WRAPPER(VECTOR_ELT)
MAGCOMP_R_WRAPPER(STRING_ELT)
MAGCOMP_R_WRAPPER(R_FlushConsole)
#undef MAGCOMP_R_WRAPPER

inline void guarded_Rprintf(const char *message)
{
  rCall([&]() { Rprintf("%s", message); });
}

template<class... Args> void guarded_Rprintf(const char *format, Args... args)
{
  rCall([&]() { Rprintf(format, args...); });
}

} // namespace magcomp
#endif
