#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"
#include <unordered_set>
#include <iomanip> 

//int main(int ac, char** av)
int main()
{
  //if (ac != 2) error_macro("expecting 1 argument: out_lgrngn parent dir")

  std::string
    //dir = string(av[1]), 
    h5  = "out_blk_1m";

  blitz::firstIndex i;
  blitz::secondIndex j;
  blitz::Range all = blitz::Range::all();
  auto n = h5n(h5);

  for (int at = 0; at < n["t"]; ++at) // TODO: mark what time does it actually mean!
  {
    for (auto &plt : std::unordered_set<std::string>({"th", "rv", "rc", "rr", "w", "u"}))
    {
      std::cout << at * n["outfreq"] << " : " << plt << std::endl;
      Gnuplot gp;
      init(gp, h5 + ".plot/" + plt + "/" + zeropad(at * n["outfreq"]) + ".png", 1, 1, n); 

      if (plt == "rc")
      {
	std::string title = "cloud water mixing ratio [g/kg]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = h5load(h5, "rc", at * n["outfreq"]);
        gp << "set cbrange [0:1.4]\n";
        plot(gp, tmp * 1000.);
      }
      else if (plt == "rr")
      {
	gp << "set logscale cb\n";
	std::string title = "rain water mixing ratio [g/kg]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = h5load(h5, "rr", at * n["outfreq"]);
        plot(gp, tmp * 1000.);
	gp << "unset logscale cb\n";
      }

      else if (plt == "rv")
      {   
	std::string title = "water vapour mixing ratio [g/kg]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = h5load(h5, "rv", at * n["outfreq"]);
        gp << "set cbrange [1:14]\n";
        plot(gp, tmp * 1000.);
      }   

      else if (plt == "th")
      {   
	std::string title = "dry air potential temperature [K]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = h5load(h5, "th", at * n["outfreq"]);
        gp << "set cbrange [288:306]\n";
        plot(gp, tmp);
      }   

      else if (plt == "u")
      {   
	std::string title = "velocity in x [m/s]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = h5load(h5, "u", at * n["outfreq"]);
        plot(gp, tmp);
      }   

      else if (plt == "w")
      {   
	std::string title = "velocity in z [m/s]"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = h5load(h5, "w", at * n["outfreq"]);
        plot(gp, tmp);
      }   

      else if (plt == "RH")
      {   
	std::string title = "relative humidity"; 
	gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
        auto tmp = h5load(h5, "RH", at * n["outfreq"]);
        plot(gp, tmp);
      }   
    } // var loop
  } // time loop
} // main
