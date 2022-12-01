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
      H = 0,
      qqZZ = 1,
      ggZZ = 2,
      rare = 3,
      ZX = 4,
      WZ = 5,
      TotalMC = 6,
      Data = 7,   
      MAX_NUM_OF_PROCESSES
   };
   
   enum _flavour
	{
		ele = 0,
		mu = 1,
		tauE = 2,
		tauMu = 3,
		tauTau = 4,
		MAX_NUM_OF_FLAVOURS
		
	};
	
   enum _final_state
	{
		fs4mu = 0,
		fs2muetau = 1,
		fs2mumutau = 2,
		fs2mutautau = 3,
		fs2e2mu = 4,
		fs2eetau = 5,
		fs2emutau = 6,
		fs2etautau = 7,
		fs4l = 8,
		MAX_NUM_OF_FINAL_STATES
	};
   
   enum _ZL_final_state
   {
       fs2mue = 0,
       fs2mumu = 1,
       fs2mutau = 2,
       fs2ee = 3,
       fs2emu = 4,
       fs2etau = 5,
       fs3l = 6,
       MAX_NUM_OF_ZL_FINAL_STATES
   };

   enum _gen_final_state
   {
        gfs4mu = 0,
        gfs2mu2tau = 1,
        gfs2e2mu = 2,
        gfs2e2tau = 3,
        gfsbkg = 4,
        gfsall = 5,
        MAX_NUM_OF_GEN_FINAL_STATES
   };

   enum _eta_bins
	{
		EB = 0,
		EE = 1,
		MAX_NUM_OF_ETA_BINS
	};
   
   enum _regions_OS
	{
		reg2P2F = 0,
		reg3P1F = 1,
		regOS   = 2,
		MAX_NUM_OF_REGIONS_OS
	};
   
   enum _regions_SS
	{
		regZLL = 0,
		MAX_NUM_OF_REGIONS_SS
	};
   
   enum _fake_rates
	{
		corrected = 0 ,
		uncorrected = 1,
		MAX_NUM_OF_FAKE_RATES
		
	};
	
	enum _FR_variations
	{
		nominal = 0,
		Up      = 1,
		Dn      = 2,
		MAX_NUM_OF_FR_VARIATIONS
	};
	
	enum _Z_mass_windows
	{
		_Window5_10 = 0,//_40_MZ1_120 = 0 ,
		_Window2_5  = 1,//_MZ1mMZtrue_7 = 1,
		_Window0_2  = 2,//_60_MZ1_120 = 2,
		//_Window0_2  = 3,//_MZ1EmMZtrue_5 = 3,
		MAX_NUM_OF_Z_MASS_WINDOWS
	};
   static const int num_of_processes         = MAX_NUM_OF_PROCESSES;
   static const int num_of_flavours          = MAX_NUM_OF_FLAVOURS;
   static const int num_of_final_states      = MAX_NUM_OF_FINAL_STATES;
   static const int num_of_zl_final_states   = MAX_NUM_OF_ZL_FINAL_STATES;
   static const int num_of_gen_final_states  = MAX_NUM_OF_GEN_FINAL_STATES;
   static const int num_of_eta_bins          = MAX_NUM_OF_ETA_BINS;
   static const int num_of_regions_os        = MAX_NUM_OF_REGIONS_OS;
   static const int num_of_regions_ss        = MAX_NUM_OF_REGIONS_SS;
   static const int num_of_fake_rates        = MAX_NUM_OF_FAKE_RATES;
   static const int num_of_fr_variations     = MAX_NUM_OF_FR_VARIATIONS;
   static const int num_of_z_mass_windows    = MAX_NUM_OF_Z_MASS_WINDOWS;

   private:
      
};
#endif
