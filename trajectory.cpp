/*  Trajectory-Calculator
 *
 *  Files: Trajector.cpp
 *
 *  Professor: P. Atrey
 *
 *  Author: Robert Bisewski
 *  Date: May 5, 2012
 */

//////////////////////
//INCLUDED LIBRARIES//
//////////////////////
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>


////////////////////////
//NAMESPACE DEFINITION//
////////////////////////
using namespace std;


//////////////////////////////////////
//GLOBAL VARIABLES ARE DECLARED HERE//
//////////////////////////////////////

const unsigned int fontsize  = 5;  //font size of labels
const unsigned int fontstyle = 7;  //font style of labels

const double alpha =  0.55;  //weighted value
const double pi    =  3.14;  //the value of pi
const double mu    =  0.50;  //the value of mu
const double sigma =  0.12;  //the value of sigma (0.12)
      double H_0   =  1.60;  //initial height
const double g     =  9.81;  //acceleration due to gravity
const double R_m   = 40.00;  //max LOS of camera
const double fps   = 50.00;  //camera refresh rate

const double V_R[] = {0.00, 1.25, 2.00, 2.50};                        //actor velocities
const double D_c[] = {5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0}; //distance coordinate

vector<double> temp; //temporarily stores data
string comparing;    //stores the comparing factor


///////////////////////////////////////////////////////////////////////////
//GLOBAL STRUCTURES: structures used globally in throughout this program.//
///////////////////////////////////////////////////////////////////////////

//stores weapon systems
class weapon
{
public:

	string name;         //weapon name
	double V_0;          //initial velocity
	double D_E;          //effective range
	double pm;           //projectile mass
	double E_0;          //initial energy
	double max;          //max range limitation (where zero = no limit enforced)
	double leth;         //lethality (as a percentage)
	double feas;         //feasibility
	vector<double> weff; //specific weapon efficiencies vector

	//constructor for weapon class
	weapon(string n, double v, double d, double m, double mx, double l, double f) {
		name = n;
		V_0  = v;
		D_E  = d;
		pm   = m;
		max  = mx;
		leth = l;
		feas = f;
		E_0  = (m * v * v) / 2;
	}

	//divides the efficiencies present in 'weff' by another given vector
	void divideBy(vector<double> vd) {
		for (unsigned int a = 0; a < weff.size(); a++) {
			if (vd.at(a) != 0) weff.at(a) = vd.at(a) / weff.at(a);
			else weff.at(a) = 0;
		}
    }
};


//////////////////////////////////////////////////////////////////////////
//GLOBAL FUNCTIONS: functions that are used to determine the efficiency.//
//////////////////////////////////////////////////////////////////////////

//calculates the height factor H_f
double height_factor(double D_c, double V_0)
{
	double hfz = 1 - ( (g * D_c * D_c) / (2 * V_0 * V_0 * H_0) );
	if (hfz > 1) hfz = 1;
	if (hfz < 0) hfz = 0;
	return hfz;
}

//calculates the velocity factor V_f
double velocity_factor(double D_c, double V_0)
{
	double vfz = 1 - ( (g * D_c) / (V_0 * V_0) );
	if (vfz > 1) vfz = 1;
	if (vfz < 0) vfz = 0;
	return vfz;
}

//calculates the range factor R_f
double range_factor(double D_E)
{
	double rfz = (D_E >= R_m) ? 1 : D_E / R_m;
	return rfz;
}

//calculates the time factor T_f
double time_factor(double V_R)
{
	double vr2 = V_R * V_R;
	double tfz = (fps - vr2) / fps;

	if (V_R < 1) tfz = 1.0;
	if (fps <= vr2) tfz = 0.0;

	return tfz;
}

//calculates the energy factor E_f
double energy_factor(double V_0, double pm, double E_0, double X)
{
	double V_x = V_0 - ((g * X) / V_0);
	double efy = (pm * V_x * V_x) / 2;
	double efz = 1 - ((E_0 - efy) / E_0);

	return efz;
}

//Gaussian probability density function
double gaussian_pdf(double damage)
{
	double gx = 0.00; //stores the exponent value
	double gy = 0.00; //stores the gaussian value

	//calculates the gaussian
	gx = (-1 * pow((damage - mu), 2)) / (2 * sigma);
	gy = (1 / (sqrt(23 * pi * sigma))) * exp(gx);

	//checks that the result is between 0 and 1
	if (gy > 1) gy = 1;
	if (gy < 0) gy = 0;

	//returns end result
	return gy;
}

//Safe exponent function that always returns a number between 0 and 1
double exp_safe(double damage)
{
    double gx = 0.00; //stores the exponent value
	double gy = 0.00; //stores the gaussian value

	//calculates the gaussian
	gx = (-1 * (damage - mu) * (damage - mu)) / (2 * sigma);
	gy = exp(gx);

	//checks that the result is between 0 and 1
	if (gy > 1) gy = 1;
	if (gy < 0) gy = 0;

	//returns end result
	return gy;
}

//DEPRECATED: determines the efficiencies of a weapon for a specific scenario
vector<double> efficiency_calc(weapon w, double V)
{
	double H_f, V_f, R_f, T_f, E_f; //factors
	double E_s, eff;                //variables used to store efficiencies
	vector<double> E_vec;           //effeciency vectors

	//calculates range and time factors
	R_f = range_factor(w.D_E);
	T_f = time_factor(V);

	//determines height/velocity factors, then the efficiencies
	for (int i = 7; i >= 0; i--) {
		H_f = height_factor(D_c[i], w.V_0);
		V_f = velocity_factor(D_c[i], w.V_0);
		E_f = energy_factor(w.V_0, w.pm, w.E_0, D_c[i]);
		E_s = (H_f * R_f * V_f * T_f);

		//calculates efficiencies
		if (E_s < 0) E_s = 0;
		eff = (alpha * E_s) + ((1 - alpha) * E_f);

		//compenstates efficiencies based on lethality
		eff *= (1 - ((w.leth-1) * 0.1));

		//adds calculated value to vector
		E_vec.push_back(eff);
	}

	//returns the completed vector
	return E_vec;
}

//determines the utility of a weapon for a specific scenario
vector<double> utility_calc(weapon w, double V, string factor, double fvalue)
{
	double H_f, R_f, T_f, E_f; //factors
	double D, U;               //variables used to damage & utility
	vector<double> U_vec;      //effeciency vectors

	//factor weights:   H_f      R_f     T_f     E_f     L_f
    const double wt[] = {0.600,  0.100,  0.100,  0.100,  0.100};

	//calculates range and time factors
	R_f = range_factor(w.D_E);
	T_f = time_factor(V);

	//determines height/velocity factors, then the utility
	for (int i = 0; i <= 7; i++) {
		//calculates H_f, E_f
		H_f = height_factor(D_c[i], w.V_0);
		E_f = energy_factor(w.V_0, w.pm, w.E_0, D_c[i]);

		//determines percentage, then calculates damage
		if (factor.compare("H") == 0) H_f = fvalue;
		if (factor.compare("R") == 0) R_f = fvalue;
		if (factor.compare("T") == 0) T_f = fvalue;
		if (factor.compare("E") == 0) E_f = fvalue;
		D = (H_f * wt[0]) + (R_f * wt[1]) + (T_f * wt[2]) + (E_f * wt[3]) + (w.leth * wt[4]);

		//using the damage, determines f(F,D)
		U = w.feas * exp_safe(D);

		//alters utility based on maximum range limitation
		if ((w.max > 0) && (w.max < D_c[i])) U = 0.0;

		//adds calculated damage value to vector
		U_vec.push_back(U);
	}

	//returns the completed vector
	return U_vec;
}

//prints the utility result table for a specific scenario
void print_table(vector<weapon> listW, int scenario, string factor)
{
	unsigned int i,j,k; //for-loop counters

	//determines if the given scenario exists
	if ((scenario < 0) || (scenario >= 4)) {

		cout << "Error: Scenario #" << scenario + 1 << " does not exist!" << endl << endl;

	} else {

		//checks if factor override is in effect
		if (factor.compare("") == 0) {

			//prints the scenario number and actor velocity
			cout << "Scenario #" << scenario + 1 << ", V_R = " << V_R[scenario] << endl;

			//prints the weapon names
			for (i = 0; i < listW.size(); i++) {
				cout << listW.at(i).name << ",";
				temp = utility_calc(listW.at(i), V_R[scenario],"",0);
				listW.at(i).weff = temp;
				temp.clear();
			}

			//determines efficiencies for each of the weapons, then prints them
			cout << endl;
			for (i = 0; i < 8; i++) {
				for (j = 0; j < listW.size(); j++) {
					cout << setprecision(4) << listW.at(j).weff.at(i) << ",";
				}
				cout << endl;
			}
			cout << endl;

		} else {

			for (k = 1; k < 5; k++) {
				//prints the scenario number and actor velocity
				cout << "Scenario #" << scenario + 1 << ", V_R = " << V_R[scenario] << ", " << factor << "_f = " << k * 0.25 << endl;

				//prints the weapon names
				for (i = 0; i < listW.size(); i++) {
					cout << listW.at(i).name << ",";
					temp = utility_calc(listW.at(i), V_R[scenario],factor,k*0.25);
					listW.at(i).weff = temp;
					temp.clear();
				}

				//determines efficiencies for each of the weapons, then prints them
				cout << endl;
				for (i = 0; i < 8; i++) {
					for (j = 0; j < listW.size(); j++) {
						cout << setprecision(4) << listW.at(j).weff.at(i) << ",";
					}
					cout << endl;
				}
				cout << endl;
			}
		}
	}
}

//prints the results as a Scilab .sc file
void print_scilab(vector<weapon> listW, int scenario, string factor)
{
	unsigned int i,j,k;                                                   //for-loop counters
	string x_label = "X";           /*"Distance (m)"*/                    //X Label
	string y_label = "U";           /*"  Utility   "*/                    //Y Label
	string L[] = {"b>-","b<-","b*-","bd-","bx-","b.-","bs-","b^-","b--"}; //line sharp array

	//determines if the given scenario exists
	if ((scenario < 0) || (scenario >= 4)) {

		cout << "Error: Scenario #" << scenario + 1 << " does not exist!" << endl << endl;

	} else {

		//checks if factor override is in effect
		if (factor.compare("") == 0) {

			//determines the weapon damages
			for (i = 0; i < listW.size(); i++) {
				temp = utility_calc(listW.at(i), V_R[scenario],"",0);
				listW.at(i).weff = temp;
				temp.clear();
			}

			//prints the value of the y-axis
			cout << "y" << scenario + 1 << "={5,10,15,20,25,30,35,40};" << endl;

			//determines efficiencies for each of the weapons, then prints them
			for (j = 0; j < listW.size(); j++) {
				cout << "x" << scenario + 1 << j << "={";
				for (i = 0; i < 8; i++) {
					cout << setprecision(4) << listW.at(j).weff.at(i);
					if (i != 7) cout << ",";
				}
				cout << "};" << endl;
			}
			cout << "plot(";
			for (j = 0; j < listW.size(); j++) {
				cout << "y" << scenario + 1 << ",x" << scenario + 1 << j << ",\"" << L[j] << "\"";
				if (j != listW.size() - 1) cout << ",";
			}
			cout << ");" << endl << "a=gca();" << endl;
			cout << "a.data_bounds = [5,0;40,1.0];" << endl;
			cout << "a.font_size=" << fontsize << ";" << endl;

			//prints the scenario number and actor velocity as per the title
			cout << "a.x_label.font_size="  << fontsize  << ";" << endl;
			cout << "a.x_label.font_style=" << fontstyle << ";" << endl;
			cout << "a.y_label.font_size="  << fontsize  << ";" << endl;
            cout << "a.y_label.font_style=" << fontstyle << ";" << endl;
			cout << "a.x_label.text='" << x_label << "';" << endl;
			cout << "a.y_label.text='" << y_label << "';" << endl;
			cout << "title('Scenario #" << scenario + 1 << ", V_R = " << V_R[scenario] << "','fontsize'," << fontsize << ");" << endl;

			//prints the graph legend
			/*cout << "legend(";
			for (i = 0; i < listW.size(); i++) {cout << "\"" << listW.at(i).name << "\""; if (i != listW.size() - 1) cout << ",";}
			cout << ",-1);" << endl << endl;*/

		} else if ((factor.compare("") != 0) && (factor.compare("Iso") != 0)) {

			for (k = 1; k < 5; k++) {

				//creates the subplot segments
				cout << "subplot(22" << k << ");" << endl;

				//determines the weapon damages
				for (i = 0; i < listW.size(); i++) {
					temp = utility_calc(listW.at(i), V_R[scenario],factor,k*0.25);
					listW.at(i).weff = temp;
					temp.clear();
				}

				//prints the value of the y-axis
				cout << "y" << scenario + 1 << "={5,10,15,20,25,30,35,40};" << endl;

				//determines efficiencies for each of the weapons, then prints them
				for (j = 0; j < listW.size(); j++) {
					cout << "x" << scenario + 1 << j << "={";
					for (i = 0; i < 8; i++) {
						cout << setprecision(4) << listW.at(j).weff.at(i);
						if (i != 7) cout << ",";
					}
					cout << "};" << endl;
				}
				cout << "plot(";
				for (j = 0; j < listW.size(); j++) {
					cout << "y" << scenario + 1 << ",x" << scenario + 1 << j << ",\"" << L[j] << "\"";
					if (j != listW.size() - 1) cout << ",";
				}
				cout << ");" << endl << "a=gca();" << endl;
				cout << "a.data_bounds = [5,0;40,1.0];" << endl;
				cout << "a.font_size=" << fontsize << ";" << endl;

				//prints the scenario number and actor velocity as per the title
			    cout << "a.x_label.font_size="  << fontsize  << ";" << endl;
                cout << "a.x_label.font_style=" << fontstyle << ";" << endl;
			    cout << "a.y_label.font_size="  << fontsize  << ";" << endl;
                cout << "a.y_label.font_style=" << fontstyle << ";" << endl;
				cout << "a.x_label.text='" << x_label << "';" << endl;
				cout << "a.y_label.text='" << y_label << "';" << endl;
				cout << "title('Scenario #" << scenario + 1 << ", V_R = " << V_R[scenario] << ", " << factor << "_f = " << k * 0.25 << "','fontsize'," << fontsize << ");" << endl;
			}
		} else {
			for (k = 1; k < 5; k++) {

				//creates the subplot segments
				cout << "subplot(22" << k << ");" << endl;

				//divides the individual factor by the overall utility
				for (i = 0; i < listW.size(); i++) {
					//calculates the initial damage
					temp = utility_calc(listW.at(i), V_R[scenario],"",0);
					listW.at(i).weff = temp;
			        temp.clear();

					//determines the final isolation effect
					temp = utility_calc(listW.at(i), V_R[scenario],comparing,k*0.25);
					listW.at(i).divideBy(temp);
					temp.clear();
				}

				//prints the value of the y-axis
				cout << "y" << scenario + 1 << "={5,10,15,20,25,30,35,40};" << endl;

				//determines efficiencies for each of the weapons, then prints them
				for (j = 0; j < listW.size(); j++) {
					cout << "x" << scenario + 1 << j << "={";
					for (i = 0; i < 8; i++) {
						cout << setprecision(4) << listW.at(j).weff.at(i);
						if (i != 7) cout << ",";
					}
					cout << "};" << endl;
				}
				cout << "plot(";
				for (j = 0; j < listW.size(); j++) {
					cout << "y" << scenario + 1 << ",x" << scenario + 1 << j << ",\"" << L[j] << "\"";
					if (j != listW.size() - 1) cout << ",";
				}
				cout << ");" << endl << "a=gca();" << endl;
				cout << "a.data_bounds = [5,0;40,3.0];" << endl;
				cout << "a.font_size=" << fontsize << ";" << endl;

				//prints the scenario number and actor velocity as per the title
			    cout << "a.x_label.font_size="  << fontsize  << ";" << endl;
			    cout << "a.x_label.font_style=" << fontstyle << ";" << endl;
			    cout << "a.y_label.font_size="  << fontsize  << ";" << endl;
                cout << "a.y_label.font_style=" << fontstyle << ";" << endl;
				cout << "a.x_label.text='" << x_label << "';" << endl;
				cout << "a.y_label.text='" << y_label << "';" << endl;
				cout << "title('Scenario #" << scenario + 1 << ", V_R = " << V_R[scenario] << ", " << factor << " = " << k * 0.25;
                cout << ", Comparing: " << comparing << "_f','fontsize'," << fontsize << ");" << endl;
			}
		}

	}

}

//prints the individual factors versus X
void print_scilab_factor(vector<weapon> listW, int scenario, string factor)
{
	unsigned int i,j;                                                     //for-loop counters
	string x_label = "X";           /*"Distance (m)"*/                    //X Label
	string y_label = factor;        /*"   Factor   "*/                    //Y Label
	string L[] = {"b>-","b<-","b*-","bd-","bx-","b.-","bs-","b^-","b--"}; //line sharp array

	//determines if the given scenario exists
	if ((scenario < 0) || (scenario >= 4)) {

		cout << "Error: Scenario #" << scenario + 1 << " does not exist!" << endl << endl;

	} else {

		//prints the value of the y-axis
		cout << "y" << scenario + 1 << "={5,10,15,20,25,30,35,40};" << endl;

		//determines the corresponding values for each factor
		for (j = 0; j < listW.size(); j++) {
			cout << "x" << scenario + 1 << j << "={";
			for (i = 0; i < 8; i++) {
				weapon w = listW.at(j);
				cout << setprecision(4);
				if (factor.compare("H") == 0) cout << height_factor(D_c[i],w.V_0);
		        if (factor.compare("R") == 0) cout << range_factor(w.D_E);
		        if (factor.compare("T") == 0) cout << time_factor(V_R[scenario]);
				if (factor.compare("E") == 0) cout << energy_factor(w.V_0,w.pm,w.E_0,D_c[i]);
				if (i != 7) cout << ",";
			}
			cout << "};" << endl;
		}

		//prepares the data for plotting
		cout << "plot(";
		for (j = 0; j < listW.size(); j++) {
			cout << "y" << scenario + 1 << ",x" << scenario + 1 << j << ",\"" << L[j] << "\"";
			if (j != listW.size() - 1) cout << ",";
		}
		cout << ");" << endl << "a=gca();" << endl;
		cout << "a.data_bounds = [5,0;40,1.0];" << endl;
		cout << "a.font_size=" << fontsize << ";" << endl;

		//prints the scenario number and actor velocity as per the title
        cout << "a.x_label.font_size="  << fontsize  << ";" << endl;
        cout << "a.x_label.font_style=" << fontstyle << ";" << endl;
        cout << "a.y_label.font_size="  << fontsize  << ";" << endl;
        cout << "a.y_label.font_style=" << fontstyle << ";" << endl;
		cout << "a.x_label.text='" << x_label << "';" << endl;
		cout << "a.y_label.text='" << y_label << "_f';" << endl;
		cout << "title('Scenario #" << scenario + 1 << ", V_R = " << V_R[scenario] << "','fontsize'," << fontsize << ");" << endl;
	}
}

//prints the comparison of factor vs utility
void print_scilab_fact_vs_util(vector<weapon> listW, int scenario, string factor, unsigned int x)
{
	unsigned int i,j;                                                     //for-loop counters
	string x_label = factor;         /*"   Utility  "*/                   //X Label
	string y_label = "U";            /*"   Factor   "*/                   //Y Label
	double Y[] = {0.2,0.4,0.6,0.8,1.0};                                   //Y Values
	string L[] = {"b>-","b<-","b*-","bd-","bx-","b.-","bs-","b^-","b--"}; //line sharp array

	//determines if the given scenario exists
	if ((scenario < 0) || (scenario >= 4)) {

		cout << "Error: Scenario #" << scenario + 1 << " does not exist!" << endl << endl;

	} else {

		//prints the value of the y-axis
		cout << "y" << scenario + 1 << "={0.2,0.4,0.6,0.8,1.0};" << endl;

		//determines the corresponding values for each factor
		for (j = 0; j < listW.size(); j++) {
			cout << "x" << scenario + 1 << j << "={";
			for (i = 0; i < 5; i++) {
				temp = utility_calc(listW.at(j), V_R[scenario],factor,Y[i]);
				cout << setprecision(4) << temp[x];
				temp.clear();
				if (i != 7) cout << ",";
			}
			cout << "};" << endl;
		}

		//prepares the data for plotting
		cout << "plot(";
		for (j = 0; j < listW.size(); j++) {
			cout << "y" << scenario + 1 << ",x" << scenario + 1 << j << ",\"" << L[j] << "\"";
			if (j != listW.size() - 1) cout << ",";
		}
		cout << ");" << endl << "a=gca();" << endl;
		cout << "a.data_bounds = [0.2,0;1,1];" << endl;
		cout << "a.font_size=" << fontsize << ";" << endl;

		//prints the scenario number and actor velocity as per the title
        cout << "a.x_label.font_size="  << fontsize  << ";" << endl;
        cout << "a.x_label.font_style=" << fontstyle << ";" << endl;
        cout << "a.y_label.font_size="  << fontsize  << ";" << endl;
        cout << "a.y_label.font_style=" << fontstyle << ";" << endl;
		cout << "a.x_label.text='" << x_label << "_f';" << endl;
		cout << "a.y_label.text='" << y_label << "';" << endl;
		cout << "title('Scenario #" << scenario + 1 << ", V_R = " << V_R[scenario] << ", X = " << D_c[x] << "','fontsize'," << fontsize << ");" << endl;
	}
}

//prints the comparison of height vs utility
void print_scilab_height_vs_util(vector<weapon> listW, int scenario, string factor, double fvalue, unsigned int x)
{
	unsigned int i,j;                                                     //for-loop counters
	string x_label = "Height";       /*"   Height    "*/                  //X Label
	string y_label = "U";            /*"   Utility   "*/                  //Y Label
	double Y[] = {0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2};                       //Y Values
	string L[] = {"b>-","b<-","b*-","bd-","bx-","b.-","bs-","b^-","b--"}; //line sharp array

	//determines if the given scenario exists
	if ((scenario < 0) || (scenario >= 4)) {

		cout << "Error: Scenario #" << scenario + 1 << " does not exist!" << endl << endl;

	} else {

		//prints the value of the y-axis
		cout << "y" << scenario + 1 << "={0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2};" << endl;

		//determines the corresponding values for each factor
		for (j = 0; j < listW.size(); j++) {
			cout << "x" << scenario + 1 << j << "={";
			for (i = 0; i < 8; i++) {
                H_0 = Y[i];
				temp = utility_calc(listW.at(j), V_R[scenario],factor,fvalue);
				cout << setprecision(4) << temp[x];
				temp.clear();
				cout << ",";
			}
			cout << "};" << endl;
		}

		//prepares the data for plotting
		cout << "plot(";
		for (j = 0; j < listW.size(); j++) {
			cout << "y" << scenario + 1 << ",x" << scenario + 1 << j << ",\"" << L[j] << "\"";
			if (j != listW.size() - 1) cout << ",";
		}
		cout << ");" << endl << "a=gca();" << endl;
		cout << "a.data_bounds = [0.4,0;3.2,1];" << endl;
		cout << "a.font_size=" << fontsize << ";" << endl;

		//prints the scenario number and actor velocity as per the title
        cout << "a.x_label.font_size="  << fontsize  << ";" << endl;
        cout << "a.x_label.font_style=" << fontstyle << ";" << endl;
        cout << "a.y_label.font_size="  << fontsize  << ";" << endl;
        cout << "a.y_label.font_style=" << fontstyle << ";" << endl;
        ///cout << "a.x_ticks=tlist(['ticks' 'locations' 'labels'],(0:8)',['0.0';'0.4';'0.8';'1.2';'1.6';'2.0';'2.4';'2.8';'3.2']);" << endl;
		cout << "a.x_label.text='" << x_label << "';" << endl;
		cout << "a.y_label.text='" << y_label << "';" << endl;
		cout << "title('Scenario #" << scenario + 1 << ", V_R = " << V_R[scenario] << ", X = " << D_c[x] << ", " << factor << "_f = " << fvalue << "','fontsize'," << fontsize << ");" << endl;
	}

    //resets the H_0 back to original value
	H_0 = 1.6;
}

/////////////////////////////////////////////////
//MAIN: The 'main' of this program begins here.//
/////////////////////////////////////////////////
int main () {

	//variable declaration
	vector<weapon> list_of_weapons; //stores all of the weapons to be used

	//defines attributes:   weapon name,  speed,  range,  weight,  max,  lethal, feas
	weapon Water_cannon( "Water Cannon",   70.0,   28.0,     1.0,    0,   0.01,  0.66);
	weapon Pepper_spray( "Pepper Spray",   70.0,   28.0,     1.0,    0,   0.01,  0.66);
	weapon Netgun(       "Netgun",         60.0,   18.0,  1000.0,    0,   0.02,  0.16);
	weapon Tear_gas(     "Tear Gas",       70.0,   28.0,     1.0,    0,   0.15,  0.66);
	weapon Crossbow(     "Crossbow",       74.0,    9.1,    32.0,    0,   0.16,  0.35);
	weapon Riot_gun(     "Riot Gun",       76.0,   27.0,   135.0,    0,   0.19,  0.21);
	weapon Taser(        "Taser",          55.0,    5.5,     1.6,   11,   0.20,  1.00);
	weapon Revolver(     "Revolver",      240.0,   82.0,    10.2,    0,   0.25,  1.00);
	weapon Sar(          "SA-Rifle",      850.0,  457.0,    34.0,    0,   0.35,  0.27);

	//adds the define weapons to the "list_of_weapons" vector
	list_of_weapons.push_back(Water_cannon);
	list_of_weapons.push_back(Netgun);
	list_of_weapons.push_back(Pepper_spray);
	list_of_weapons.push_back(Tear_gas);
	list_of_weapons.push_back(Riot_gun);
	list_of_weapons.push_back(Taser);
	list_of_weapons.push_back(Crossbow);
	list_of_weapons.push_back(Revolver);
	list_of_weapons.push_back(Sar);

	//displays the utility result tables for the various scenarios
	cout << "subplot(221)" << endl;
	print_scilab(list_of_weapons,0,"");
	cout << "subplot(222)" << endl;
	print_scilab(list_of_weapons,1,"");
	cout << "subplot(223)" << endl;
	print_scilab(list_of_weapons,2,"");
	cout << "subplot(224)" << endl;
	print_scilab(list_of_weapons,3,"");

    //scripts for the H,R,T,E factors
    print_scilab(list_of_weapons,0,"H");
    print_scilab(list_of_weapons,0,"R");
    print_scilab(list_of_weapons,0,"T");
    print_scilab(list_of_weapons,0,"E");

	//prints the scilab script for the isolated scenario
	comparing = "E";
	print_scilab(list_of_weapons,0,"Iso");

	//compares each factor to X
	cout << "subplot(221)" << endl;
	print_scilab_factor(list_of_weapons,0,"H");
	cout << "subplot(222)" << endl;
	print_scilab_factor(list_of_weapons,0,"R");
	cout << "subplot(223)" << endl;
	print_scilab_factor(list_of_weapons,0,"T");
	cout << "subplot(224)" << endl;
	print_scilab_factor(list_of_weapons,0,"E");

	//attempts to compare a factor to the overal utility
    cout << "subplot(221)" << endl;
	print_scilab_fact_vs_util(list_of_weapons,0,"H",5);
	cout << "subplot(222)" << endl;
	print_scilab_fact_vs_util(list_of_weapons,0,"R",5);
	cout << "subplot(223)" << endl;
	print_scilab_fact_vs_util(list_of_weapons,0,"T",5);
	cout << "subplot(224)" << endl;
	print_scilab_fact_vs_util(list_of_weapons,0,"E",5);

	//plots height versus utility
	cout << "subplot(221)" << endl;
	print_scilab_height_vs_util(list_of_weapons, 0, "E", 0.25, 5);
	cout << "subplot(222)" << endl;
	print_scilab_height_vs_util(list_of_weapons, 0, "E", 0.50, 5);
	cout << "subplot(223)" << endl;
	print_scilab_height_vs_util(list_of_weapons, 0, "E", 0.75, 5);
	cout << "subplot(224)" << endl;
	print_scilab_height_vs_util(list_of_weapons, 0, "E", 1.00, 5);

	//this ends the program
	return 0;
}
