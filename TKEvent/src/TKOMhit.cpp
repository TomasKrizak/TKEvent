// TK headers
#include "TKOMhit.h"

ClassImp(TKOMhit);

void TKOMhit::set_SWCR_xyz()
{
	//mainwall IT
	if(OM_num < 260) 
	{
		SWCR[0] = 0;
		SWCR[1] = -1;
		SWCR[2] = OM_num / 13;
		SWCR[3] = OM_num % 13;
	}
	//mainwall FR
	else if(OM_num < 520)
	{
		SWCR[0] = 1;
		SWCR[1] = -1;
		SWCR[2] = (OM_num - 260) / 13;
		SWCR[3] = (OM_num - 260) % 13;
	}
	//Xcalo IT	
	else if(OM_num < 584)
	{
		SWCR[0] = 0;
		SWCR[1] = (OM_num - 520) / 32;
		SWCR[2] = ((OM_num - 520) / 16) % 2;
		SWCR[3] = (OM_num -520) % 16;
	}
	//Xcalo FR
	else if(OM_num < 648)
	{
		SWCR[0] = 1;
		SWCR[1] = (OM_num - 520 - 64) / 32;
		SWCR[2] = ((OM_num - 520 - 64) / 16) % 2;
		SWCR[3] = (OM_num -520 - 64) % 16;
	}
	//GVeto IT
	else if(OM_num < 680)
	{
		SWCR[0] = 0;
		SWCR[1] = (OM_num - 520 - 128) / 16;
		SWCR[2] = (OM_num - 520 - 128) % 16;
		SWCR[3] = -1;
	}
	//GVeto FR
	else if(OM_num < 712)
	{
		SWCR[0] = 1;
		SWCR[1] = (OM_num - 520 - 128 - 32) / 16;
		SWCR[2] = (OM_num - 520 - 128 - 32) % 16;
		SWCR[3] = -1;
	}
	
	int OM_type;
	
	if(OM_num < 520)
	{
		OM_type = 1302;
	}
	else if(OM_num < 648)
	{
		OM_type = 1232;
	}
	else
	{
		OM_type = 1252;
	}

	switch(OM_type)
	{
		case 1302: //MW
			if(SWCR[0] == 1)
				xyz[0] = 532.0;
			else
				xyz[0] = -532.0;
				xyz[1] = ((double)SWCR[2]- 9.5) * 259.0;
				xyz[2] = ((double)SWCR[3] - 6) * 259.0;
				
			break;
			
		case 1232: //XW
			if(SWCR[1] == 1)
				xyz[1] = 2580.5;
			else
				xyz[1] = -2580.5;
				
			if(SWCR[0] == 1)
			{
				if(SWCR[2] == 1)
					xyz[0] = 333.0;
				else
					xyz[0] = 130.0;
			}
			else
			{
				if(SWCR[2] == 1)
					xyz[0] = -333.0;
				else
					xyz[0] = -130.0;
			}
			
			xyz[2] = ((double)SWCR[3] - 7.5) * 212.0;
			
			break;
			
		case 1252: //GV
			if(SWCR[0] == 1)
				xyz[0] = 213.5;
			else
				xyz[0] = -213.5;
			if(SWCR[1] == 1)
				xyz[2] = 1625.0;
			else
				xyz[2] = -1625.0;
			if(SWCR[2] > 7)
				xyz[1] = 161.0 + (((double)SWCR[2]-8) * 311.5);
			else
				xyz[1] = -161.0 + (((double)SWCR[2]-7) * 311.5);
			break;	
// MIRO: PRIDAT DEFAULT PRE PRIPAD CHYBNEHO VSTUPU?
	}
}

TKOMhit::TKOMhit()
{
}

TKOMhit::TKOMhit(int _OM_num)
{
// MIRO: DAŤ PODMIENKU ČO KONTROLUJE ČÍSLO OM!!
//	if(_cell_num > 2033 || _cell_num < 0) std::cout << "ERROR TKtrhit::TKtrhit(int): " << _cell_num << " is not valid tracker cell number" << std::endl; 

	OM_num = _OM_num;
	set_SWCR_xyz();
}

TKOMhit::TKOMhit(int _SWCR[4])
{
	// auto detect MW
	if (_SWCR[0] != -1 && 
	    _SWCR[1] == -1 && 
	    _SWCR[2] != -1 && 
	    _SWCR[3] != -1 )
	{
		OM_num = 260*_SWCR[0] + 13*_SWCR[2] + _SWCR[3];
	}
	// auto detect XW
	else if (_SWCR[0] != -1 && 
		 _SWCR[1] != -1 && 
		 _SWCR[2] != -1 && 
		 _SWCR[3] != -1 )
	{
		OM_num = 520 + 64*_SWCR[0] + 32*_SWCR[1] + 16*_SWCR[2] + _SWCR[3];
	}
	// auto detect GV
	else if (_SWCR[0] != -1 && 
		 _SWCR[1] != -1 && 
		 _SWCR[2] != -1 && 
		 _SWCR[3] == -1 )
	{
		OM_num = 520 + 128 + 32*_SWCR[0] + 16*_SWCR[1] + _SWCR[2];
	}
	else 
	{
// MIRO: ŠTANDARDIZOVAŤ VÝPIS ČO KONTROLUJE ČÍSLO OM!!
		std::cout << "warning: " << _SWCR[0] << "." << _SWCR[1] << "." << _SWCR[2] << "." << _SWCR[3] << " is not valid OM" << std::endl;
	}
	
// MIRO: PREHODNOTIT ROZDELENIE SWCR A XYZ
	set_SWCR_xyz();
}

TKOMhit::TKOMhit(int _OM_num, bool _HT, int64_t _OM_TDC, int16_t _OM_pcell)
{
// MIRO: DAŤ PODMIENKU ČO KONTROLUJE ČÍSLO OM!!
//	if(_cell_num > 2033 || _cell_num < 0) std::cout << "ERROR TKtrhit::TKtrhit(int): " << _cell_num << " is not valid tracker cell number" << std::endl; 

	OM_num = _OM_num;
	set_SWCR_xyz();
	
	set_HT(_HT);
	set_OM_TDC(_OM_TDC);
	set_OM_pcell(_OM_pcell);
}

TKOMhit::TKOMhit(int _SWCR[4], bool _HT, int64_t _OM_TDC, int16_t _OM_pcell)
{
	// auto detect MW
	if (_SWCR[0] != -1 && 
	    _SWCR[1] == -1 && 
	    _SWCR[2] != -1 && 
	    _SWCR[3] != -1 )
	{
		OM_num = 260*_SWCR[0] + 13*_SWCR[2] + _SWCR[3];
	}
	// auto detect XW
	else if (_SWCR[0] != -1 && 
		 _SWCR[1] != -1 && 
		 _SWCR[2] != -1 && 
		 _SWCR[3] != -1 )
	{
		OM_num = 520 + 64*_SWCR[0] + 32*_SWCR[1] + 16*_SWCR[2] + _SWCR[3];
	}
	// auto detect GV
	else if (_SWCR[0] != -1 && 
		 _SWCR[1] != -1 && 
		 _SWCR[2] != -1 && 
		 _SWCR[3] == -1 )
	{
		OM_num = 520 + 128 + 32*_SWCR[0] + 16*_SWCR[1] + _SWCR[2];
	}
	else 
	{
// MIRO: ŠTANDARDIZOVAŤ VÝPIS ČO KONTROLUJE ČÍSLO OM!!
		std::cout << "warning: " << _SWCR[0] << "." << _SWCR[1] << "." << _SWCR[2] << "." << _SWCR[3] << " is not valid OM" << std::endl;
	}
	
// MIRO: PREHODNOTIT ROZDELENIE SWCR A XYZ
	set_SWCR_xyz();
	
	set_HT(_HT);
	set_OM_TDC(_OM_TDC);
	set_OM_pcell(_OM_pcell);
}

TKOMhit::~TKOMhit()
{
}

void TKOMhit::set_OM_num(int _OM_num)
{
	OM_num   = _OM_num;
}

void TKOMhit::set_HT(bool _HT)
{
	HT    = _HT;
}

void TKOMhit::set_OM_TDC(int64_t _OM_TDC)
{
	OM_TDC   = _OM_TDC;
}

void TKOMhit::set_OM_pcell(int16_t _OM_pcell)
{
	OM_pcell = _OM_pcell;
}

int TKOMhit::get_OM_num()
{
	return OM_num;
}

int TKOMhit::get_SWCR(char _SWCR_n)
{	
	switch (_SWCR_n)
	{
		case 's':	return SWCR[0];
			  	break;
		case 'S':	return SWCR[0];
			  	break;
		case 'w':	return SWCR[1];
			  	break;
		case 'W':	return SWCR[1];
			  	break;
		case 'c':	return SWCR[2];
			  	break;
		case 'C':	return SWCR[2];
			  	break;
		case 'r':	return SWCR[3];
			  	break;
		case 'R':	return SWCR[3];
			  	break;
		default:	std::cout << "ERROR in int TKOMhit::get_SWCR(char): " << _SWCR_n << " is not a valid argument value! " << std::endl;  
				return NULL;
	}
}

double TKOMhit::get_xyz(char _xyz_n)
{	
	switch (_xyz_n)
	{
		case 'x':	return xyz[0];
			  	break;
		case 'X':	return xyz[0];
			  	break;
		case 'y':	return xyz[1];
			  	break;
		case 'Y':	return xyz[1];
			  	break;
		case 'z':	return xyz[2];
			  	break;
		case 'Z':	return xyz[2];
			  	break;
		default:	std::cout << "ERROR in double TKOMhit::get_xyz(char): " << _xyz_n << " is not a valid argument value! " << std::endl;  
				return NULL;
	}
}

bool    TKOMhit::is_HT()
{
	return HT;
}

int64_t TKOMhit::get_OM_TDC()
{
	return OM_TDC;
}

int64_t TKOMhit::get_OM_pcell()
{
	return OM_pcell;
}

void TKOMhit::print()
{
	std::cout << "	OM "       << OM_num 
	     	  << ", is HT: " << HT 
	     	  << ", tdc: "   << OM_TDC
	     	  << ", pcell: " << OM_pcell << std::endl;
}
