// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdOdIlibdITKOMhitdict

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
#include "TKOMhit.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TKOMhit(void *p = 0);
   static void *newArray_TKOMhit(Long_t size, void *p);
   static void delete_TKOMhit(void *p);
   static void deleteArray_TKOMhit(void *p);
   static void destruct_TKOMhit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TKOMhit*)
   {
      ::TKOMhit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TKOMhit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TKOMhit", ::TKOMhit::Class_Version(), "TKOMhit.h", 10,
                  typeid(::TKOMhit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TKOMhit::Dictionary, isa_proxy, 4,
                  sizeof(::TKOMhit) );
      instance.SetNew(&new_TKOMhit);
      instance.SetNewArray(&newArray_TKOMhit);
      instance.SetDelete(&delete_TKOMhit);
      instance.SetDeleteArray(&deleteArray_TKOMhit);
      instance.SetDestructor(&destruct_TKOMhit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TKOMhit*)
   {
      return GenerateInitInstanceLocal((::TKOMhit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TKOMhit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TKOMhit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TKOMhit::Class_Name()
{
   return "TKOMhit";
}

//______________________________________________________________________________
const char *TKOMhit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKOMhit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TKOMhit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKOMhit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TKOMhit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKOMhit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TKOMhit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKOMhit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TKOMhit::Streamer(TBuffer &R__b)
{
   // Stream an object of class TKOMhit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TKOMhit::Class(),this);
   } else {
      R__b.WriteClassBuffer(TKOMhit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TKOMhit(void *p) {
      return  p ? new(p) ::TKOMhit : new ::TKOMhit;
   }
   static void *newArray_TKOMhit(Long_t nElements, void *p) {
      return p ? new(p) ::TKOMhit[nElements] : new ::TKOMhit[nElements];
   }
   // Wrapper around operator delete
   static void delete_TKOMhit(void *p) {
      delete ((::TKOMhit*)p);
   }
   static void deleteArray_TKOMhit(void *p) {
      delete [] ((::TKOMhit*)p);
   }
   static void destruct_TKOMhit(void *p) {
      typedef ::TKOMhit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TKOMhit

namespace {
  void TriggerDictionaryInitialization_TKOMhitdict_Impl() {
    static const char* headers[] = {
"TKOMhit.h",
0
    };
    static const char* includePaths[] = {
"/sps/nemo/sw/BxCppDev/opt/root-6.16.00/include/root",
"/pbs/home/t/tkrizak/supernemo-commissioning/TKEvent/TKEvent/include/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TKOMhitdict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$TKOMhit.h")))  TKOMhit;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TKOMhitdict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TKOMhit.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TKOMhit", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TKOMhitdict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TKOMhitdict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TKOMhitdict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TKOMhitdict() {
  TriggerDictionaryInitialization_TKOMhitdict_Impl();
}
