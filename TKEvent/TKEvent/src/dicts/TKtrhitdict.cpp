// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdOdIlibdITKtrhitdict
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
#include "TKtrhit.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_TKtrhit(void *p = 0);
   static void *newArray_TKtrhit(Long_t size, void *p);
   static void delete_TKtrhit(void *p);
   static void deleteArray_TKtrhit(void *p);
   static void destruct_TKtrhit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TKtrhit*)
   {
      ::TKtrhit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TKtrhit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TKtrhit", ::TKtrhit::Class_Version(), "TKtrhit.h", 11,
                  typeid(::TKtrhit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TKtrhit::Dictionary, isa_proxy, 4,
                  sizeof(::TKtrhit) );
      instance.SetNew(&new_TKtrhit);
      instance.SetNewArray(&newArray_TKtrhit);
      instance.SetDelete(&delete_TKtrhit);
      instance.SetDeleteArray(&deleteArray_TKtrhit);
      instance.SetDestructor(&destruct_TKtrhit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TKtrhit*)
   {
      return GenerateInitInstanceLocal((::TKtrhit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TKtrhit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TKtrhit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TKtrhit::Class_Name()
{
   return "TKtrhit";
}

//______________________________________________________________________________
const char *TKtrhit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKtrhit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TKtrhit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKtrhit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TKtrhit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKtrhit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TKtrhit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKtrhit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TKtrhit::Streamer(TBuffer &R__b)
{
   // Stream an object of class TKtrhit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TKtrhit::Class(),this);
   } else {
      R__b.WriteClassBuffer(TKtrhit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TKtrhit(void *p) {
      return  p ? new(p) ::TKtrhit : new ::TKtrhit;
   }
   static void *newArray_TKtrhit(Long_t nElements, void *p) {
      return p ? new(p) ::TKtrhit[nElements] : new ::TKtrhit[nElements];
   }
   // Wrapper around operator delete
   static void delete_TKtrhit(void *p) {
      delete ((::TKtrhit*)p);
   }
   static void deleteArray_TKtrhit(void *p) {
      delete [] ((::TKtrhit*)p);
   }
   static void destruct_TKtrhit(void *p) {
      typedef ::TKtrhit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TKtrhit

namespace {
  void TriggerDictionaryInitialization_TKtrhitdict_Impl() {
    static const char* headers[] = {
"TKtrhit.h",
0
    };
    static const char* includePaths[] = {
"/home/tomas/Programs/root/install/include/",
"/home/tomas/TKEvent/TKEvent/include/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TKtrhitdict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$TKtrhit.h")))  TKtrhit;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TKtrhitdict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "TKtrhit.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"TKtrhit", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TKtrhitdict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TKtrhitdict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TKtrhitdict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TKtrhitdict() {
  TriggerDictionaryInitialization_TKtrhitdict_Impl();
}
