// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdOdIlibdITKtimerdict
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
#include "TKtimer.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_TKtimer(void *p = 0);
   static void *newArray_TKtimer(Long_t size, void *p);
   static void delete_TKtimer(void *p);
   static void deleteArray_TKtimer(void *p);
   static void destruct_TKtimer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TKtimer*)
   {
      ::TKtimer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TKtimer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TKtimer", ::TKtimer::Class_Version(), "TKtimer.h", 11,
                  typeid(::TKtimer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TKtimer::Dictionary, isa_proxy, 4,
                  sizeof(::TKtimer) );
      instance.SetNew(&new_TKtimer);
      instance.SetNewArray(&newArray_TKtimer);
      instance.SetDelete(&delete_TKtimer);
      instance.SetDeleteArray(&deleteArray_TKtimer);
      instance.SetDestructor(&destruct_TKtimer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TKtimer*)
   {
      return GenerateInitInstanceLocal((::TKtimer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TKtimer*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TKtimer::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TKtimer::Class_Name()
{
   return "TKtimer";
}

//______________________________________________________________________________
const char *TKtimer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKtimer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TKtimer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKtimer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TKtimer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKtimer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TKtimer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKtimer*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TKtimer::Streamer(TBuffer &R__b)
{
   // Stream an object of class TKtimer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TKtimer::Class(),this);
   } else {
      R__b.WriteClassBuffer(TKtimer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TKtimer(void *p) {
      return  p ? new(p) ::TKtimer : new ::TKtimer;
   }
   static void *newArray_TKtimer(Long_t nElements, void *p) {
      return p ? new(p) ::TKtimer[nElements] : new ::TKtimer[nElements];
   }
   // Wrapper around operator delete
   static void delete_TKtimer(void *p) {
      delete ((::TKtimer*)p);
   }
   static void deleteArray_TKtimer(void *p) {
      delete [] ((::TKtimer*)p);
   }
   static void destruct_TKtimer(void *p) {
      typedef ::TKtimer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TKtimer

namespace {
  void TriggerDictionaryInitialization_TKtimerdict_Impl() {
    static const char* headers[] = {
"TKtimer.h",
0
    };
    static const char* includePaths[] = {
"/home/tomas/Programs/root/install/include/",
"/home/tomas/TKEvent/TKEvent/include/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TKtimerdict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$TKtimer.h")))  TKtimer;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TKtimerdict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "TKtimer.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"TKtimer", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TKtimerdict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TKtimerdict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TKtimerdict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TKtimerdict() {
  TriggerDictionaryInitialization_TKtimerdict_Impl();
}
