#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"
#include <unordered_set>
#include <iomanip> 

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: dir with data")

  //std::string h5 = string(av[1]);
  for (auto &aer : std::unordered_set<std::string>({"080"}))//, "110", "140", "170", "200", "230"}))
  {
    //std::string h5 = "out_blk_2m_piggy_" + aer;
    std::string h5 = "out_blk_2m_driver_080_rf01";

    blitz::firstIndex i;
    blitz::secondIndex j;
    blitz::Range all = blitz::Range::all();

    auto n = h5n(h5);

    for (int at = 0; at < n["t"]; ++at) // TODO: mark what time does it actually mean!
    {
      //for (auto &plt : std::unordered_set<std::string>({"th", "rv", "rc", "rr", "u", "w", "one", "thousand"}))
      for (auto &plt : std::unordered_set<std::string>({"th", "rv", "rc", "rr", "nr", "nc", "u", "w", "one", "thousand"}))
      {
        std::cout << at * n["outfreq"] << " : " << plt << std::endl;
        Gnuplot gp;
        init(gp, h5 + ".plot/" + plt + "/" + zeropad(at * n["outfreq"]) + ".png", 1, 1, n); 

        if (plt == "one")
        {
          std::string title = "one [-]"; 
          gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
          auto tmp = h5load(h5, "one", at * n["outfreq"]);
          //gp << "set cbrange [-2:2]\n";
          plot(gp, tmp);
        }
        if (plt == "thousand")
        {
          std::string title = "thousand [-]"; 
          gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
          auto tmp = h5load(h5, "thousand", at * n["outfreq"]);
          gp << "set cbrange [999.8:1000.25]\n";
          plot(gp, tmp);
        }

        if (plt == "rc")
        {
          std::string title = "cloud water mixing ratio [g/kg]"; 
          gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
          auto tmp = h5load(h5, "rc", at * n["outfreq"]);
          gp << "set cbrange [0:3.]\n";
          plot(gp, tmp * 1000.);
        }

        else if (plt == "rr")
        {
          std::string title = "rain water mixing ratio [g/kg]"; 
          gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
          auto tmp = h5load(h5, "rr", at * n["outfreq"]);
          gp << "set logscale cb\n";
          gp << "set cbrange [0.001:10.]\n";
          plot(gp, tmp * 1000.);
          gp << "unset logscale cb\n";
        }

        if (plt == "nc")
        {
          std::string title = "cloud droplet spec conc [#/mg]"; 
          gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
          auto tmp = h5load(h5, "nc", at * n["outfreq"]);
          gp << "set cbrange [0:250]\n";
          plot(gp, tmp * 1e-6);
        }
        if (plt == "nr")
        {
          std::string title = "rain drop spec conc [#/mg]"; 
          gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
          auto tmp = h5load(h5, "nr", at * n["outfreq"]);
          gp << "set cbrange [0.001:10]\n";
          gp << "set logscale cb\n";
          plot(gp, tmp * 1e-6);
          gp << "unset logscale cb\n";
        }
        else if (plt == "rv")
        {   
          std::string title = "water vapour mixing ratio [g/kg]"; 
          gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
          auto tmp = h5load(h5, "rv", at * n["outfreq"]);
          gp << "set cbrange [2:16]\n";
          plot(gp, tmp * 1000.);
        }   

        else if (plt == "th")
        {   
          std::string title = "dry air potential temperature [K]"; 
          gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
          auto tmp = h5load(h5, "th", at * n["outfreq"]);
          gp << "set cbrange [288:310]\n";
          plot(gp, tmp);
        }   

        else if (plt == "u")
        {   
          std::string title = "velocity in x [m/s]"; 
          gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
          auto tmp = h5load(h5, "u", at * n["outfreq"]);
          //gp << "set cbrange [-0.01:0.01]\n";
          plot(gp, tmp);
        }   

        else if (plt == "w")
        {   
          std::string title = "velocity in z [m/s]"; 
          gp << "set title '" + title + " t = " << std::fixed << std::setprecision(2) << (double(at) * n["outfreq"] * n["dt"] / 60.) << "min'\n";
          auto tmp = h5load(h5, "w", at * n["outfreq"]);
          //gp << "set cbrange [-0.008:0.008]\n";
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
  } //aer loop
} // main