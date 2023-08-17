#ifndef Settings_h
#define Settings_h

using namespace std;

class Settings
{

public:
	
	Settings();
	~Settings();
   
    enum _process
    {
      Htt = 0,
	  Hmm = 1,
	  Hww = 2,
	  diBoson = 3,
	  triBoson = 4,
	  top = 5,
	  fake = 6,
	  EWKZee = 6,
	  EWKZtt = 7,
	  EWKZmm = 8,
	  DYee = 9,
	  DYtt = 10,
	  DYmm = 11,
      TotalMC = 12,
      Data = 13,   
      MAX_NUM_OF_PROCESSES
    };
	
    enum _final_state
	{
		fsmumu = 0,
		fsetau = 1,
		fsmutau = 2,
		fstautau = 3,
		fsemu = 4,
		fs2l = 5,
		MAX_NUM_OF_FINAL_STATES
	};
	
	enum _category
	{
		catGGH = 0,
		catVBF = 1,
		catBoost = 2,
		catAll = 3,
		MAX_NUM_OF_CATEGORIES
	};

	enum _FR_variations
	{
		nominal = 0,
		Up      = 1,
		Dn      = 2,
		MAX_NUM_OF_FR_VARIATIONS
	};
	
   static const int num_of_processes         = MAX_NUM_OF_PROCESSES;
   static const int num_of_final_states      = MAX_NUM_OF_FINAL_STATES;
   static const int num_of_categories		 = MAX_NUM_OF_CATEGORIES;
   static const int num_of_fr_variations     = MAX_NUM_OF_FR_VARIATIONS;

   private:
      
};
#endif
