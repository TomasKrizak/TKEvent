// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdOdIlibdITKtrackdict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "TKtrack.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TKtrack(void *p = 0);
   static void *newArray_TKtrack(Long_t size, void *p);
   static void delete_TKtrack(void *p);
   static void deleteArray_TKtrack(void *p);
   static void destruct_TKtrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TKtrack*)
   {
      ::TKtrack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TKtrack >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TKtrack", ::TKtrack::Class_Version(), "TKtrack.h", 14,
                  typeid(::TKtrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TKtrack::Dictionary, isa_proxy, 4,
                  sizeof(::TKtrack) );
      instance.SetNew(&new_TKtrack);
      instance.SetNewArray(&newArray_TKtrack);
      instance.SetDelete(&delete_TKtrack);
      instance.SetDeleteArray(&deleteArray_TKtrack);
      instance.SetDestructor(&destruct_TKtrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TKtrack*)
   {
      return GenerateInitInstanceLocal((::TKtrack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TKtrack*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TKtrack::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TKtrack::Class_Name()
{
   return "TKtrack";
}

//______________________________________________________________________________
const char *TKtrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKtrack*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TKtrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKtrack*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TKtrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKtrack*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TKtrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKtrack*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TKtrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class TKtrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TKtrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(TKtrack::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TKtrack(void *p) {
      return  p ? new(p) ::TKtrack : new ::TKtrack;
   }
   static void *newArray_TKtrack(Long_t nElements, void *p) {
      return p ? new(p) ::TKtrack[nElements] : new ::TKtrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_TKtrack(void *p) {
      delete ((::TKtrack*)p);
   }
   static void deleteArray_TKtrack(void *p) {
      delete [] ((::TKtrack*)p);
   }
   static void destruct_TKtrack(void *p) {
      typedef ::TKtrack current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TKtrack

namespace {
  void TriggerDictionaryInitialization_TKtrackdict_Impl() {
    static const char* headers[] = {
"TKtrack.h",
0
    };
    static const char* includePaths[] = {
"/sps/nemo/sw/BxCppDev/opt/root-6.16.00/include/root",
"/pbs/home/t/tkrizak/supernemo-commissioning/TKEvent/TKEvent/include/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TKtrackdict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$TKtrack.h")))  TKtrack;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TKtrackdict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TKtrack.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TKtrack", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TKtrackdict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TKtrackdict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TKtrackdict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TKtrackdict() {
  TriggerDictionaryInitialization_TKtrackdict_Impl();
}
