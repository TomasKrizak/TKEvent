// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdOdIlibdITKEventdict
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
#include "TKEvent.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_TKEvent(void *p = 0);
   static void *newArray_TKEvent(Long_t size, void *p);
   static void delete_TKEvent(void *p);
   static void deleteArray_TKEvent(void *p);
   static void destruct_TKEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TKEvent*)
   {
      ::TKEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TKEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TKEvent", ::TKEvent::Class_Version(), "TKEvent.h", 49,
                  typeid(::TKEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TKEvent::Dictionary, isa_proxy, 4,
                  sizeof(::TKEvent) );
      instance.SetNew(&new_TKEvent);
      instance.SetNewArray(&newArray_TKEvent);
      instance.SetDelete(&delete_TKEvent);
      instance.SetDeleteArray(&deleteArray_TKEvent);
      instance.SetDestructor(&destruct_TKEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TKEvent*)
   {
      return GenerateInitInstanceLocal((::TKEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TKEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TKEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TKEvent::Class_Name()
{
   return "TKEvent";
}

//______________________________________________________________________________
const char *TKEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TKEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TKEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TKEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TKEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TKEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TKEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class TKEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TKEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(TKEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TKEvent(void *p) {
      return  p ? new(p) ::TKEvent : new ::TKEvent;
   }
   static void *newArray_TKEvent(Long_t nElements, void *p) {
      return p ? new(p) ::TKEvent[nElements] : new ::TKEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_TKEvent(void *p) {
      delete ((::TKEvent*)p);
   }
   static void deleteArray_TKEvent(void *p) {
      delete [] ((::TKEvent*)p);
   }
   static void destruct_TKEvent(void *p) {
      typedef ::TKEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TKEvent

namespace ROOT {
   static TClass *vectorlETKtrhitmUgR_Dictionary();
   static void vectorlETKtrhitmUgR_TClassManip(TClass*);
   static void *new_vectorlETKtrhitmUgR(void *p = 0);
   static void *newArray_vectorlETKtrhitmUgR(Long_t size, void *p);
   static void delete_vectorlETKtrhitmUgR(void *p);
   static void deleteArray_vectorlETKtrhitmUgR(void *p);
   static void destruct_vectorlETKtrhitmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TKtrhit*>*)
   {
      vector<TKtrhit*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TKtrhit*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TKtrhit*>", -2, "vector", 386,
                  typeid(vector<TKtrhit*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETKtrhitmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TKtrhit*>) );
      instance.SetNew(&new_vectorlETKtrhitmUgR);
      instance.SetNewArray(&newArray_vectorlETKtrhitmUgR);
      instance.SetDelete(&delete_vectorlETKtrhitmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETKtrhitmUgR);
      instance.SetDestructor(&destruct_vectorlETKtrhitmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TKtrhit*> >()));

      ::ROOT::AddClassAlternate("vector<TKtrhit*>","std::vector<TKtrhit*, std::allocator<TKtrhit*> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TKtrhit*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETKtrhitmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TKtrhit*>*)0x0)->GetClass();
      vectorlETKtrhitmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETKtrhitmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETKtrhitmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TKtrhit*> : new vector<TKtrhit*>;
   }
   static void *newArray_vectorlETKtrhitmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TKtrhit*>[nElements] : new vector<TKtrhit*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETKtrhitmUgR(void *p) {
      delete ((vector<TKtrhit*>*)p);
   }
   static void deleteArray_vectorlETKtrhitmUgR(void *p) {
      delete [] ((vector<TKtrhit*>*)p);
   }
   static void destruct_vectorlETKtrhitmUgR(void *p) {
      typedef vector<TKtrhit*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TKtrhit*>

namespace ROOT {
   static TClass *vectorlETKtrajectorymUgR_Dictionary();
   static void vectorlETKtrajectorymUgR_TClassManip(TClass*);
   static void *new_vectorlETKtrajectorymUgR(void *p = 0);
   static void *newArray_vectorlETKtrajectorymUgR(Long_t size, void *p);
   static void delete_vectorlETKtrajectorymUgR(void *p);
   static void deleteArray_vectorlETKtrajectorymUgR(void *p);
   static void destruct_vectorlETKtrajectorymUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TKtrajectory*>*)
   {
      vector<TKtrajectory*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TKtrajectory*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TKtrajectory*>", -2, "vector", 386,
                  typeid(vector<TKtrajectory*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETKtrajectorymUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TKtrajectory*>) );
      instance.SetNew(&new_vectorlETKtrajectorymUgR);
      instance.SetNewArray(&newArray_vectorlETKtrajectorymUgR);
      instance.SetDelete(&delete_vectorlETKtrajectorymUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETKtrajectorymUgR);
      instance.SetDestructor(&destruct_vectorlETKtrajectorymUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TKtrajectory*> >()));

      ::ROOT::AddClassAlternate("vector<TKtrajectory*>","std::vector<TKtrajectory*, std::allocator<TKtrajectory*> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TKtrajectory*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETKtrajectorymUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TKtrajectory*>*)0x0)->GetClass();
      vectorlETKtrajectorymUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETKtrajectorymUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETKtrajectorymUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TKtrajectory*> : new vector<TKtrajectory*>;
   }
   static void *newArray_vectorlETKtrajectorymUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TKtrajectory*>[nElements] : new vector<TKtrajectory*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETKtrajectorymUgR(void *p) {
      delete ((vector<TKtrajectory*>*)p);
   }
   static void deleteArray_vectorlETKtrajectorymUgR(void *p) {
      delete [] ((vector<TKtrajectory*>*)p);
   }
   static void destruct_vectorlETKtrajectorymUgR(void *p) {
      typedef vector<TKtrajectory*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TKtrajectory*>

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
   static TClass *vectorlETKclustermUgR_Dictionary();
   static void vectorlETKclustermUgR_TClassManip(TClass*);
   static void *new_vectorlETKclustermUgR(void *p = 0);
   static void *newArray_vectorlETKclustermUgR(Long_t size, void *p);
   static void delete_vectorlETKclustermUgR(void *p);
   static void deleteArray_vectorlETKclustermUgR(void *p);
   static void destruct_vectorlETKclustermUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TKcluster*>*)
   {
      vector<TKcluster*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TKcluster*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TKcluster*>", -2, "vector", 386,
                  typeid(vector<TKcluster*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETKclustermUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TKcluster*>) );
      instance.SetNew(&new_vectorlETKclustermUgR);
      instance.SetNewArray(&newArray_vectorlETKclustermUgR);
      instance.SetDelete(&delete_vectorlETKclustermUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETKclustermUgR);
      instance.SetDestructor(&destruct_vectorlETKclustermUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TKcluster*> >()));

      ::ROOT::AddClassAlternate("vector<TKcluster*>","std::vector<TKcluster*, std::allocator<TKcluster*> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TKcluster*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETKclustermUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TKcluster*>*)0x0)->GetClass();
      vectorlETKclustermUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETKclustermUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETKclustermUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TKcluster*> : new vector<TKcluster*>;
   }
   static void *newArray_vectorlETKclustermUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TKcluster*>[nElements] : new vector<TKcluster*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETKclustermUgR(void *p) {
      delete ((vector<TKcluster*>*)p);
   }
   static void deleteArray_vectorlETKclustermUgR(void *p) {
      delete [] ((vector<TKcluster*>*)p);
   }
   static void destruct_vectorlETKclustermUgR(void *p) {
      typedef vector<TKcluster*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TKcluster*>

namespace ROOT {
   static TClass *vectorlETKOMhitmUgR_Dictionary();
   static void vectorlETKOMhitmUgR_TClassManip(TClass*);
   static void *new_vectorlETKOMhitmUgR(void *p = 0);
   static void *newArray_vectorlETKOMhitmUgR(Long_t size, void *p);
   static void delete_vectorlETKOMhitmUgR(void *p);
   static void deleteArray_vectorlETKOMhitmUgR(void *p);
   static void destruct_vectorlETKOMhitmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TKOMhit*>*)
   {
      vector<TKOMhit*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TKOMhit*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TKOMhit*>", -2, "vector", 386,
                  typeid(vector<TKOMhit*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETKOMhitmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TKOMhit*>) );
      instance.SetNew(&new_vectorlETKOMhitmUgR);
      instance.SetNewArray(&newArray_vectorlETKOMhitmUgR);
      instance.SetDelete(&delete_vectorlETKOMhitmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETKOMhitmUgR);
      instance.SetDestructor(&destruct_vectorlETKOMhitmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TKOMhit*> >()));

      ::ROOT::AddClassAlternate("vector<TKOMhit*>","std::vector<TKOMhit*, std::allocator<TKOMhit*> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TKOMhit*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETKOMhitmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TKOMhit*>*)0x0)->GetClass();
      vectorlETKOMhitmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETKOMhitmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETKOMhitmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TKOMhit*> : new vector<TKOMhit*>;
   }
   static void *newArray_vectorlETKOMhitmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TKOMhit*>[nElements] : new vector<TKOMhit*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETKOMhitmUgR(void *p) {
      delete ((vector<TKOMhit*>*)p);
   }
   static void deleteArray_vectorlETKOMhitmUgR(void *p) {
      delete [] ((vector<TKOMhit*>*)p);
   }
   static void destruct_vectorlETKOMhitmUgR(void *p) {
      typedef vector<TKOMhit*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TKOMhit*>

namespace {
  void TriggerDictionaryInitialization_TKEventdict_Impl() {
    static const char* headers[] = {
"TKEvent.h",
0
    };
    static const char* includePaths[] = {
"/home/tomas/Programs/root/install/include/",
"/home/tomas/TKEvent/TKEvent/include/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TKEventdict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$TKEvent.h")))  TKEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TKEventdict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "TKEvent.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"TKEvent", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TKEventdict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TKEventdict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TKEventdict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TKEventdict() {
  TriggerDictionaryInitialization_TKEventdict_Impl();
}
