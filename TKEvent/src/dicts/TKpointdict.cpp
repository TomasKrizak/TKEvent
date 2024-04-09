// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdOdIlibdITKpointdict
#define R__NO_DEPRECATION

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

// Header files passed as explicit arguments
#include "TKpoint.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_TKpoint(void *p = 0);
   static void *newArray_TKpoint(Long_t size, void *p);
   static void delete_TKpoint(void *p);
   static void deleteArray_TKpoint(void *p);
   static void destruct_TKpoint(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TKpoint*)
   {
      ::TKpoint *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TKpoint >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TKpoint", ::TKpoint::Class_Version(), "TKpoint.h", 11,
                  typeid(::TKpoint), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TKpoint::Dictionary, isa_proxy, 4,
                  sizeof(::TKpoint) );
      instance.SetNew(&new_TKpoint);
      instance.SetNewArray(&newArray_TKpoint);
      instance.SetDelete(&delete_TKpoint);
      instance.SetDeleteArray(&deleteArray_TKpoint);
      instance.SetDestructor(&destruct_TKpoint);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TKpoint*)
   {
      return GenerateInitInstanceLocal((::TKpoint*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TKpoint*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TKpoint::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TKpoint::Class_Name()
{
   return "TKpoint";
}

//______________________________________________________________________________
const char *TKpoint::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKpoint*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TKpoint::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKpoint*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TKpoint::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKpoint*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TKpoint::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKpoint*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TKpoint::Streamer(TBuffer &R__b)
{
   // Stream an object of class TKpoint.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TKpoint::Class(),this);
   } else {
      R__b.WriteClassBuffer(TKpoint::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TKpoint(void *p) {
      return  p ? new(p) ::TKpoint : new ::TKpoint;
   }
   static void *newArray_TKpoint(Long_t nElements, void *p) {
      return p ? new(p) ::TKpoint[nElements] : new ::TKpoint[nElements];
   }
   // Wrapper around operator delete
   static void delete_TKpoint(void *p) {
      delete ((::TKpoint*)p);
   }
   static void deleteArray_TKpoint(void *p) {
      delete [] ((::TKpoint*)p);
   }
   static void destruct_TKpoint(void *p) {
      typedef ::TKpoint current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TKpoint

namespace {
  void TriggerDictionaryInitialization_TKpointdict_Impl() {
    static const char* headers[] = {
"TKpoint.h",
0
    };
    static const char* includePaths[] = {
"/home/tomas/Programs/root/install/include/",
"/home/tomas/TKEvent/TKEvent/include/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TKpointdict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$TKpoint.h")))  TKpoint;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TKpointdict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "TKpoint.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"TKpoint", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TKpointdict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TKpointdict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TKpointdict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TKpointdict() {
  TriggerDictionaryInitialization_TKpointdict_Impl();
}
