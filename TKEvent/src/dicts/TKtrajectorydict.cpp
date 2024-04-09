// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdOdIlibdITKtrajectorydict
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
#include "TKtrajectory.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_TKtrajectory(void *p = 0);
   static void *newArray_TKtrajectory(Long_t size, void *p);
   static void delete_TKtrajectory(void *p);
   static void deleteArray_TKtrajectory(void *p);
   static void destruct_TKtrajectory(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TKtrajectory*)
   {
      ::TKtrajectory *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TKtrajectory >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TKtrajectory", ::TKtrajectory::Class_Version(), "TKtrajectory.h", 17,
                  typeid(::TKtrajectory), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TKtrajectory::Dictionary, isa_proxy, 4,
                  sizeof(::TKtrajectory) );
      instance.SetNew(&new_TKtrajectory);
      instance.SetNewArray(&newArray_TKtrajectory);
      instance.SetDelete(&delete_TKtrajectory);
      instance.SetDeleteArray(&deleteArray_TKtrajectory);
      instance.SetDestructor(&destruct_TKtrajectory);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TKtrajectory*)
   {
      return GenerateInitInstanceLocal((::TKtrajectory*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TKtrajectory*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TKtrajectory::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TKtrajectory::Class_Name()
{
   return "TKtrajectory";
}

//______________________________________________________________________________
const char *TKtrajectory::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKtrajectory*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TKtrajectory::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKtrajectory*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TKtrajectory::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKtrajectory*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TKtrajectory::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKtrajectory*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TKtrajectory::Streamer(TBuffer &R__b)
{
   // Stream an object of class TKtrajectory.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TKtrajectory::Class(),this);
   } else {
      R__b.WriteClassBuffer(TKtrajectory::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TKtrajectory(void *p) {
      return  p ? new(p) ::TKtrajectory : new ::TKtrajectory;
   }
   static void *newArray_TKtrajectory(Long_t nElements, void *p) {
      return p ? new(p) ::TKtrajectory[nElements] : new ::TKtrajectory[nElements];
   }
   // Wrapper around operator delete
   static void delete_TKtrajectory(void *p) {
      delete ((::TKtrajectory*)p);
   }
   static void deleteArray_TKtrajectory(void *p) {
      delete [] ((::TKtrajectory*)p);
   }
   static void destruct_TKtrajectory(void *p) {
      typedef ::TKtrajectory current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TKtrajectory

namespace ROOT {
   static TClass *vectorlETKtrackmUgR_Dictionary();
   static void vectorlETKtrackmUgR_TClassManip(TClass*);
   static void *new_vectorlETKtrackmUgR(void *p = 0);
   static void *newArray_vectorlETKtrackmUgR(Long_t size, void *p);
   static void delete_vectorlETKtrackmUgR(void *p);
   static void deleteArray_vectorlETKtrackmUgR(void *p);
   static void destruct_vectorlETKtrackmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TKtrack*>*)
   {
      vector<TKtrack*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TKtrack*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TKtrack*>", -2, "vector", 386,
                  typeid(vector<TKtrack*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETKtrackmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TKtrack*>) );
      instance.SetNew(&new_vectorlETKtrackmUgR);
      instance.SetNewArray(&newArray_vectorlETKtrackmUgR);
      instance.SetDelete(&delete_vectorlETKtrackmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETKtrackmUgR);
      instance.SetDestructor(&destruct_vectorlETKtrackmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TKtrack*> >()));

      ::ROOT::AddClassAlternate("vector<TKtrack*>","std::vector<TKtrack*, std::allocator<TKtrack*> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TKtrack*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETKtrackmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TKtrack*>*)0x0)->GetClass();
      vectorlETKtrackmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETKtrackmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETKtrackmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TKtrack*> : new vector<TKtrack*>;
   }
   static void *newArray_vectorlETKtrackmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TKtrack*>[nElements] : new vector<TKtrack*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETKtrackmUgR(void *p) {
      delete ((vector<TKtrack*>*)p);
   }
   static void deleteArray_vectorlETKtrackmUgR(void *p) {
      delete [] ((vector<TKtrack*>*)p);
   }
   static void destruct_vectorlETKtrackmUgR(void *p) {
      typedef vector<TKtrack*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TKtrack*>

namespace ROOT {
   static TClass *vectorlETKpointmUgR_Dictionary();
   static void vectorlETKpointmUgR_TClassManip(TClass*);
   static void *new_vectorlETKpointmUgR(void *p = 0);
   static void *newArray_vectorlETKpointmUgR(Long_t size, void *p);
   static void delete_vectorlETKpointmUgR(void *p);
   static void deleteArray_vectorlETKpointmUgR(void *p);
   static void destruct_vectorlETKpointmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TKpoint*>*)
   {
      vector<TKpoint*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TKpoint*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TKpoint*>", -2, "vector", 386,
                  typeid(vector<TKpoint*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETKpointmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TKpoint*>) );
      instance.SetNew(&new_vectorlETKpointmUgR);
      instance.SetNewArray(&newArray_vectorlETKpointmUgR);
      instance.SetDelete(&delete_vectorlETKpointmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETKpointmUgR);
      instance.SetDestructor(&destruct_vectorlETKpointmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TKpoint*> >()));

      ::ROOT::AddClassAlternate("vector<TKpoint*>","std::vector<TKpoint*, std::allocator<TKpoint*> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TKpoint*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETKpointmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TKpoint*>*)0x0)->GetClass();
      vectorlETKpointmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETKpointmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETKpointmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TKpoint*> : new vector<TKpoint*>;
   }
   static void *newArray_vectorlETKpointmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TKpoint*>[nElements] : new vector<TKpoint*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETKpointmUgR(void *p) {
      delete ((vector<TKpoint*>*)p);
   }
   static void deleteArray_vectorlETKpointmUgR(void *p) {
      delete [] ((vector<TKpoint*>*)p);
   }
   static void destruct_vectorlETKpointmUgR(void *p) {
      typedef vector<TKpoint*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TKpoint*>

namespace {
  void TriggerDictionaryInitialization_TKtrajectorydict_Impl() {
    static const char* headers[] = {
"TKtrajectory.h",
0
    };
    static const char* includePaths[] = {
"/home/tomas/Programs/root/install/include/",
"/home/tomas/TKEvent/TKEvent/include/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TKtrajectorydict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$TKtrajectory.h")))  TKtrajectory;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TKtrajectorydict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "TKtrajectory.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"TKtrajectory", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TKtrajectorydict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TKtrajectorydict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TKtrajectorydict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TKtrajectorydict() {
  TriggerDictionaryInitialization_TKtrajectorydict_Impl();
}
