#pragma once

// root stuff
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TChain.h>
#include <ROOT/RDataFrame.hxx>

// other stuff
#include <filesystem>
#include <boost/format.hpp>
#include <iostream>

// my stuff
#include "plots.cpp"
#include "../plot_style.cpp"

using namespace std;
using boost::format;

// this is just for easy printing. note that they are based on the setup.json file
const map<int, string> detector_map = {{0, "SU"}, {1, "SD"}, {2, "Det1"}, {3, "Det2"}};

void print_title(string title) {
    std::cout << "\n\033[1;32m" + title + "\033[0m" << endl;
    return;
}

// filter the input to only contain elements from a given detector
template <typename T>
vector<T> detector_filter(const vector<T> *input, const vector<int> *ID, int detector) {
    int n = ID->size(); 
    const vector<int> &IDref = *ID; 
    const vector<T> &inputref = *input;
    vector<T> result(n);
    int m = 0; // current index in result
    for (int i = 0; i < n; i++) {
        if (IDref[i] == detector) {
            result[m] = inputref[i];
            m++;
        }
    }
    result.resize(m);
    return result; // remove all the unset entries
}

// a variant of the detector filter which filters FT, BT, FI, BI, ID simultaneously. I needed this multiple times in my code, so I made this specialized function for it
tuple<vector<double>, vector<double>, vector<int>, vector<int>, vector<int>> detector_filter(const vector<double> *FT, const vector<double> *BT, const vector<int> *FI, const vector<int> *BI, const vector<int> *ID, int detector) {
    // dereference input
    int n = ID->size(); 
    const vector<int> &IDref = *ID; 
    const vector<double> &FTref = *FT;
    const vector<double> &BTref = *BT;
    const vector<int> &FIref = *FI;
    const vector<int> &BIref = *BI;

    // initialize result vectors
    vector<double> ft(n);
    vector<double> bt(n);
    vector<int> fi(n);
    vector<int> bi(n);
    vector<int> id(n);

    int m = 0; // current index in results
    for (int i = 0; i < n; i++) {
        if (IDref[i] == detector) {
            ft[m] = FTref[i];
            bt[m] = BTref[i];
            fi[m] = FIref[i];
            bi[m] = BIref[i];
            id[m] = IDref[i];
            m++;
        }
    }
    // remove all unset entries and return the tuple
    ft.resize(m); bt.resize(m); fi.resize(m); bi.resize(m); id.resize(m);
    return tuple(ft, bt, fi, bi, id);
}

// append the second input to the first
template <typename T>
void append(vector<T> *v1, vector<T> *v2) {
    v1->insert(v1->end(), v2->begin(), v2->end());
    return;
}

// append the second input to the first
template <typename T>
void append(vector<T> *v1, vector<T> v2) {
    v1->insert(v1->end(), v2.begin(), v2.end());
    return;
}

// finds coincidences between the given S3 detector and the two W1 detectors
// returns a tuple for each input. note the first two generic F and T vector inputs, which allows preserving any type of information
template <typename T>
tuple<vector<T>, vector<T>, vector<int>, vector<int>, vector<int>> find_coincidences(const vector<T> *F, const vector<T> *B, const vector<int> *FI, const vector<int> *BI, const vector<int> *ID, int det) {
    // dereference the input (does not copy them - they are references to the data itself)
    const vector<T> &Fref = *F; 
    const vector<T> &Bref = *B; 
    const vector<int> &FIref = *FI; 
    const vector<int> &BIref = *BI; 
    const vector<int> &IDref = *ID; 

    // initialize result vectors
    int n = FI->size();
    vector<T> rf(n);
    vector<T> rb(n);
    vector<int> rfi(n);
    vector<int> rbi(n);
    vector<int> rid(n); // id of the W detector

    // current index counters for the result vectors
    int m = 0;

    // checks multiple-detector coincidences between three events, and adds them to the result list if found
    // note that entries are always duplicated. this is because if a multiple-detector coincidence is found, from p1 | p2 | p3
    // one of them must have been found by another detector than the other two. without loss of generality, let this be p1. 
    // then both p1 | p2 and p1 | p3 are added, since they are two distinct coincidences that can be used for analysis
    auto check_coincidences = [&] (int i1, int i2, int i3) {
        int p1 = IDref[i1]; int p2 = IDref[i2]; int p3 = IDref[i3];
        if (p1 == det) { // S3
            if (p2 == 2 || p2 == 3) { // W1
                // p2 is from W1
                rid[m] = IDref[i2]; 
                rf[m] = Fref[i2]; // F is always taken from the W1 detectors (2 || 3)
                rfi[m] = FIref[i2]; // the associated strip number

                // p1 is from S3
                rb[m] = Bref[i1]; // we always take B from the S3 detector
                rbi[m] = BIref[i1]; // the associated strip number 
                m++; // increment all indices
            }
            if (p3 == 2 || p3 == 3) { // W1
                // p3 is from W1
                rid[m] = IDref[i3];
                rf[m] = Fref[i3]; 
                rfi[m] = FIref[i3];

                // p1 is from S3
                rb[m] = Bref[i1];
                rbi[m] = BIref[i1];
                m++;
            }
        } else if (p1 == 2 || p1 == 3) { // W1
            if (p2 == det) {
                // p1 is from W1
                rid[m] = IDref[i1]; 
                rf[m] = Fref[i1]; 
                rfi[m] = FIref[i1];

                // p2 is from S3
                rb[m] = Bref[i2];
                rbi[m] = BIref[i2];
                m++;
            }
            if (p3 == det) { // S3
                // p1 is from W1
                rid[m] = IDref[i1]; 
                rf[m] = Fref[i1]; 
                rfi[m] = FIref[i1];

                // p3 is from S3
                rb[m] = Bref[i3];
                rbi[m] = BIref[i3];
                m++;
            }
        }
        // check coincidences between p2 and p3. we have already covered all cases for p1
        if (p2 == det) { // S3
            if (p3 == 2 || p3 == 3) { // W1
                // p3 is from W1
                rid[m] = IDref[i3];
                rf[m] = Fref[i3];
                rfi[m] = FIref[i3];

                // p2 is from S3
                rb[m] = Bref[i2];
                rbi[m] = BIref[i2];
                m++;
            }
        } else if (p2 == 2 || p2 == 3) { // W1
            if (p3 == det) { // S3
                // p2 is from W1
                rid[m] = IDref[i2];
                rf[m] = Fref[i2];
                rfi[m] = FIref[i2];

                // p3 is from S3
                rb[m] = Bref[i3]; 
                rbi[m] = BIref[i3];
                m++;
            }
        }
    };

    // imagining the original data as three columns a1 | a2 | a3, it is unpacked by simply stacking these columns (a1 a2 a3)
    // thus we actually have not lost any information about the coincidences, which we use here
    n = n/3; 
    for (int i = 0; i < n; i++) {
        check_coincidences(i, n+i, 2*n+i);
    }
    rf.resize(m); rb.resize(m); rfi.resize(m); rbi.resize(m); rid.resize(m);
    return tuple(rf, rb, rfi, rbi, rid);
}

// calculate dt = (bt-ft)
vector<double> calc_dt(const vector<double> *FT, const vector<double> *BT) {
    int n = FT->size();
    const vector<double> &FTref = *FT;
    const vector<double> &BTref = *BT;
    vector<double> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = BTref[i] - FTref[i];
    }
    return result;
}

// finds coincidences between the given S3 detector and the two W1 detectors
// this is simply a wrapper for when a vector property does not depend solely on either F or B (e.g. deltaE)
// returns two vectors with the values associated with the front and back strip, respectively. 
template <typename T>
pair<vector<T>, vector<T>> find_coincidences(const vector<T> *input, const vector<int> *FI, const vector<int> *BI, const vector<int> *ID, int det) {
    vector<T> output_F, output_B;
    tie(output_F, output_B, std::ignore, std::ignore, std::ignore) = find_coincidences(input, input, FI, BI, ID, det);
    return pair(output_F, output_B);
}

// removes all entries with BT == 0, and, if det == 2 or 3, also those with FT == 0. 
void remove_0(vector<double> *FT, vector<double> *BT, vector<int> *FI, vector<int> *BI, vector<int> *ID, int det) {
    // dereference the input
    vector<double> &ft = *FT;
    vector<double> &bt = *BT;
    vector<int> &fi = *FI;
    vector<int> &bi = *BI;
    vector<int> &id = *ID;

    std::function<bool(int, int, int)> check;
    if (det == 0 || det == 1) { 
        // checks if any of the three elements have bt == 0
        check = [&] (int i1, int i2, int i3) {
            if (bt[i1] == 0) return false;
            if (bt[i2] == 0) return false;
            if (bt[i3] == 0) return false;
            return true;
        };
    } else { // also checks if any element has ft == 0
        check = [&] (int i1, int i2, int i3) {
            if (bt[i1] == 0) return false;
            if (bt[i2] == 0) return false;
            if (bt[i3] == 0) return false;

            if (ft[i1] == 0) return false;
            if (ft[i2] == 0) return false;
            if (ft[i3] == 0) return false;
            return true;};
    }

    int n = ft.size(); 
    int m = n/3;

    vector<bool> filter(n, false);
    int i2, i3;
    // loop over all decay events
    for (int i1 = 0; i1 < m; i1++) {
        i2 = i1 + m; i3 = i2 + m;
        if (check(i1, i2, i3)) {
            filter[i1] = true;
            filter[i2] = true;
            filter[i3] = true;
        }
    }

    int c = 0;
    // loop over all entries
    for (int i = 0; i < n; i++) {
        if (filter[i]) {
            ft[c] = ft[i];
            bt[c] = bt[i];
            fi[c] = fi[i];
            bi[c] = bi[i];
            id[c] = id[i];
            c++;
        }
    }
    ft.resize(c); bt.resize(c); fi.resize(c); bi.resize(c); id.resize(c);
    return;
}

// a container which allows me to modify an arbitrary number of vectors simultaneously
// this is needed for imposing cuts on all data when we don't even know what all data is (it may change in the future)
// sadly I thought of this approach too late, so it is not used as well as it could've been
// there must definitely be some way of implementing both doubles and integers as the same data type. I'm not sure how, though, and I don't want to dedicate more time to figuring out how
class data_container {
    public: 
        data_container() {}

        template <typename T>
        void check_n(string name, vector<T>* v) {
            if (n == 0) {
                n = v->size();
            } else {
                if (n != v->size()) {
                    std::cout << "\033[1;31m" << "Warning: Added vector \"" << name << "\" does not have the same size as the others!" << "\033[0m" << endl; // red colour
                }
            }
        }

        // adds a dt column to the internal data. it is corrected to be non-zero for the S3 detectors
        void add_dt() {
            vector<double>* FT = get_double("FT");
            vector<double>* BT = get_double("BT");
            vector<int>* FI = get_int("FI");
            vector<int>* BI = get_int("BI");
            vector<int>* ID = get_int("ID");

            // calculate the dt vector. note that all entries from the S3's are 0
            vector<double> dt = calc_dt(FT, BT);

            // fix dt entries for the S3 detectors
            for (int det = 0; det <= 1; det++) {
                // set up the index vector
                vector<int> index(n);
                for (int i = 0; i < n; i++) {
                    index[i] = i;
                }
                // extract the back indices of each coincidence event
                tie(std::ignore, index) = find_coincidences(&index, FI, BI, ID, det);

                // extract the other necessary stuff from the coincidences, and calculate dt
                auto[ft, bt, fi, bi, id] = find_coincidences(FT, BT, FI, BI, ID, det);
                vector<double> dt_S3 = calc_dt(&ft, &bt);

                // hist_info info(plot::x_axes::standard, plot::y_axis.at(det), "Title", "dt [ns]", "Back strip id");
                // hist2D(&dt_S3, &bi, info);

                // correct each dt entry in the main dt vector
                for (int i = 0; i < index.size(); i++) {
                    dt[index[i]] = dt_S3[i];
                }
            }
            add_data("dt", dt);
        }

        // add an int vector to the container
        void add_data(string name, vector<int> v) {
            if (check_int("name")) {
                std::cout << format("A column named %1% already exists. Overwriting.") % name << endl;
            } else {
                int_names.push_back(name);
                n_int++;
            }
            int_data[name] = v;
            check_n(name, &v);
        }

        // add a double vector to the container
        void add_data(string name, vector<double> v) {
            if (check_double("name")) {
                std::cout << format("A column named %1% already exists. Overwriting.") % name << endl;
            } else {
                double_names.push_back(name);
                n_double++;
            }
            double_data[name] = v;
            check_n(name, &v);
        }

        // check if an integer column named "name" exist
        bool check_int(string name) {
            return int_data.count(name) == 1;
        }

        // check if a double column named "name" exist
        bool check_double(string name) {
            return double_data.count(name) == 1;
        }

        // return a pointer to the int data with name "name"
        vector<int>* get_int(string name) {
            if (!check_int(name)) {
                std::cout << "\033[1;31m" << "Warning: \"" << name << "\" not found." << "\033[0m" << endl;
            }
            return &int_data.at(name);
        }

        // return a pointer to the double data with name "name"
        vector<double>* get_double(string name) {
            if (!check_double(name)) {
                std::cout << "\033[1;31m" << "Warning: \"" << name << "\" not found." << "\033[0m" << endl;
            }
            return &double_data.at(name);
        }

        int get_n() {
            return n;
        }

        void dt_filter(const vector<double> a, const vector<double> b) {
            // dereference the relevant data
            add_dt();
            const vector<double> &dt = *get_double("dt");
            const vector<int> &id = *get_int("ID"); // we will need this one often, so dereference it

            int m = n/3; // iterate through whole 12C decay events instead of alpha particles
            vector<bool> filt(n, false); // only "true" events will be kept
            auto check = [&] (int i1, int i2, int i3) {
                // check if any of the three alphas are outside the bounds
                if (!(id[i1]==-1) && !(a[id[i1]] < dt[i1] && dt[i1] < b[id[i1]])) return;
                if (!(id[i2]==-1) && !(a[id[i2]] < dt[i2] && dt[i2] < b[id[i2]])) return;
                if (!(id[i3]==-1) && !(a[id[i3]] < dt[i3] && dt[i3] < b[id[i3]])) return;

                // if so, add them all to the filter
                filt[i1] = true;
                filt[i2] = true;
                filt[i3] = true;
                return;
            };

            // iterate through each 12C decay event
            int i2, i3;
            for (int i1 = 0; i1 < m; i1++) {
                i2 = i1 + m; 
                i3 = i1 + 2*m;
                check(i1, i2, i3);
            }

            filter(&filt);
            return;
        }

        // impose a dt filter on all data, treating it as whole 12C decay events. preserving determines whether the (a1 | a2 | a3) structure should be preserved or not
        void dt_filter(const vector<double> a, const vector<double> b, bool preserving) {
            // dereference the relevant data
            add_dt();
            const vector<double> &dt = *get_double("dt");
            const vector<int> &id = *get_int("ID"); // we will need this one often, so dereference it

            int m = n/3; // iterate through whole 12C decay events instead of alpha particles
            vector<bool> filt(n, false); // only "true" events will be kept

            // helper function. checks if all three indices are within the bounds of their detectors
            std::function<void(int, int, int)> check;

            if (preserving) { // the preserving check
                check = [&] (int i1, int i2, int i3) {
                    // check if any of the three alphas are outside the bounds
                    if (!(a[id[i1]] < dt[i1] && dt[i1] < b[id[i1]])) return;
                    if (!(a[id[i2]] < dt[i2] && dt[i2] < b[id[i2]])) return;
                    if (!(a[id[i3]] < dt[i3] && dt[i3] < b[id[i3]])) return;

                    // if so, add them all to the filter
                    filt[i1] = true;
                    filt[i2] = true;
                    filt[i3] = true;
                    return;
                };
            } else { // the non-preserving check
                check = [&] (int i1, int i2, int i3) {
                    // check if each individual alpha is within the bounds, and add them to the filter if so
                    if (a[id[i1]] < dt[i1] && dt[i1] < b[id[i1]]) filt[i1] = true;
                    if (a[id[i2]] < dt[i2] && dt[i2] < b[id[i2]]) filt[i2] = true;
                    if (a[id[i3]] < dt[i3] && dt[i3] < b[id[i3]]) filt[i3] = true;
                    return;
                };
            }

            // iterate through each 12C decay event
            int i2, i3;
            for (int i1 = 0; i1 < m; i1++) {
                i2 = i1 + m; 
                i3 = i1 + 2*m;
                check(i1, i2, i3);
            }

            filter(&filt);
            return;
        }

        // impose a dt filter on all data individually. this means that we lose the (a1 | a2 | a3) structure
        // note that since this breaks the structure, the S3 histograms cannot be reconstructed after this method is called
        void dt_filter_soft(const vector<double> a, const vector<double> b) {
            const vector<double>* FT = get_double("FT");
            const vector<double>* BT = get_double("BT");
            const vector<int>* FI = get_int("FI");
            const vector<int>* BI = get_int("BI");
            const vector<int>* ID = get_int("ID");
            const vector<int> &IDref = *ID;
            vector<bool> filter(n, false);

            // calculate the dt vector. note that all entries from the S3's are 0
            vector<double> dt = calc_dt(FT, BT);

            // fix dt entries for the S3 detectors
            for (int det = 0; det <= 1; det++) {
                vector<double> ft, bt, dt_S3;
                vector<int> index(n);

                // set up the index vector
                for (int i = 0; i < n; i++) {
                    index[i] = i;
                }
                // extract the back indices of each coincidence event
                tie(std::ignore, index) = find_coincidences(&index, FI, BI, ID, det);

                tie(ft, bt, std::ignore, std::ignore, std::ignore) = find_coincidences(FT, BT, FI, BI, ID, det);
                dt_S3 = calc_dt(&ft, &bt);
                for (int i = 0; i < index.size(); i++) {
                    dt[index[i]] = dt_S3[i];
                }
            }

            int det;
            for (int i = 0; i < n; i++) {
                det = IDref[i]; // the detector which observed this event
                if (a[det] < dt[i] && dt[i] < b[det]) { // filter to events satisfying the dt cut
                    filter[i] = true;
                }
            }

            // we perform the actual filtering in-place, so we avoid having to double the memory usage
            int m; // current index in result
            for (int i = 0; i < n_int; i++) { // loop over all names in int_names
                vector<int> &col = int_data[int_names[i]];
                m = 0;
                for (int j = 0; j < n; j++) { // loop over all entries in the given vector
                    if (filter[j]) {
                        col[m] = col[j];
                        m++;
                    }
                }
                col.resize(m); // remove all unset entries
            }

            for (int i = 0; i < n_double; i++) { // loop over all vector<int> data sets
                vector<double> &col = double_data[double_names[i]];
                m = 0;
                for (int j = 0; j < n; j++) { // loop over all entries in one vector
                    if (filter[j]) {
                        col[m] = col[j];
                        m++;
                    }
                }
                col.resize(m); // remove all unset entries
            }

            n = m; // update the size of the vectors
            return;
        }

        data_container clone() {
            data_container result;
            for (int i = 0; i < n_int; i++) { // clone int data
                result.add_data(int_names[i], int_data[int_names[i]]);
            }
            for (int i = 0; i < n_double; i++) { // clone double data
                result.add_data(double_names[i], double_data[double_names[i]]);
            }
            return result;
        }

        // impose a filter on all data. note that filter must be stacked three times on top of itself
        void filter(vector<bool>* f) {
            for (int i = 0; i < n_int; i++) {
                vector<int>* col = &int_data[int_names[i]];
                filter(col, f);
            }
            for (int i = 0; i < n_double; i++) {
                vector<double>* col = &double_data[double_names[i]];
                filter(col, f);
            }
            
            // update n
            vector<bool>& filt = *f;
            n = 0;
            for (int i = 0; i < filt.size(); i++) {
                if (filt[i]) {
                    n++;
                }
            }
            
            return;
        }

    private:
        // sets of names, so we have something to iterate over
        vector<string> int_names;
        vector<string> double_names;

        // the actual data storage
        map<string, vector<int>> int_data;
        map<string, vector<double>> double_data;

        int n = 0; // size of each vector
        int n_int = 0; // number of int vectors
        int n_double = 0; // number of double vectors

        // filters a given vector
        template <typename T>
        void filter(vector<T>* vec, vector<bool>* filter) {
            vector<bool> &f = *filter;
            vector<T> &v = *vec;
            int c = 0; // current index in result
            for (int i = 0; i < v.size(); i++) {
                if (f[i]) {
                    v[c] = v[i];
                    c++;
                }
            }
            v.resize(c);
            return;
        }

};

// applies the offsets from all detectors to bt
void apply_offset(data_container* data, const vector<vector<double>> offset) {
    data->add_dt();
    const vector<int> &bi = *data->get_int("BI");
    const vector<int> &id = *data->get_int("ID");
    vector<double> &bt = *data->get_double("BT");
    
    for (int i = 0; i < bt.size(); i++) {
        bt[i] -= offset[id[i]][bi[i]-1];
    }
    return;
}

// subtracts the offsets specified in offset from each element of the vector DT
void apply_offset(vector<double> *DT, const vector<int> *BI, const vector<double> offset) {
    const vector<int> &BIref = *BI;
    vector<double> &DTref = *DT;
    for (int i = 0; i < DTref.size(); i++) {
        DTref[i] -= offset[BIref[i]-1];
    }
    return;
}

// applies deltaF and deltaB to FT and BT according to the fitted equation
// we cannot simply extract the results from the fitting class, since FT and BT would then be scrambled compared to the others
void apply_fit(vector<double> *FT, vector<double> *BT, const vector<int> *FI, const vector<int> *BI, const vector<double> doubleF, const vector<double> doubleB) {
    vector<double> &FTref = *FT;
    vector<double> &BTref = *BT;
    const vector<int> &FIref = *FI;
    const vector<int> &BIref = *BI; 

    for (int i = 0; i < FTref.size(); i++) {
        FTref[i] += doubleF[FIref[i]-1];
        BTref[i] += doubleB[BIref[i]-1];
    }
    return;
}

// applies deltaF to FT according to the fitted equation
// this is used to apply a W1 fit to the S3 data
void apply_fit(vector<double> *FT, const vector<int> *FI, const vector<int> *ID, const vector<vector<double>> deltaF) {
    vector<double> &FTref = *FT;
    const vector<int> &FIref = *FI;
    const vector<int> &IDref = *ID; 

    for (int i = 0; i < FTref.size(); i++) {
        FTref[i] += deltaF[IDref[i]][FIref[i]-1];
    }
    return;
}

// applies deltaF and deltaB to FT and BT according to the fitted equation
// this method applies the fit on data from all detectors simultaneously
void apply_fit(data_container* data, const vector<vector<double>> deltaF, const vector<vector<double>> deltaB) {
    vector<double> &FTref = *data->get_double("FT");
    vector<double> &BTref = *data->get_double("BT");
    const vector<int> &FIref = *data->get_int("FI");
    const vector<int> &BIref = *data->get_int("BI");
    const vector<int> &IDref = *data->get_int("ID");

    for (int i = 0; i < FTref.size(); i++) {
        FTref[i] += deltaF[IDref[i]][FIref[i]-1];
        BTref[i] += deltaB[IDref[i]][BIref[i]-1];
    }
    return;
}

// calculate the mean, used as the offset in this code
vector<double> calc_mean(const vector<int> *BI, const vector<double> *DT, int strips) {
    const vector<int> &BIref = *BI;
    const vector<double> &DTref = *DT;
    vector<double> sum(strips);
    vector<int> count(strips);
    
    // sum over each of the strips
    for (int i = 0; i < BIref.size(); i++) {
        int bi = BIref[i]-1;
        sum[bi] += DTref[i];
        count[bi]++;
    }
    
    // divide each sum by its number of elements (i.e. calculate the mean)
    for (int i = 0; i < strips; i++) {
        if (count[i] != 0) {
            sum[i] /= count[i];
        } else {
            sum[i] = 0;
        }
    }
    return sum;
}

// calculate the mode for each back strip
vector<double> calc_mode(const vector<int> *BI, const vector<double> *DT, int strips) {
    const vector<int>& bi = *BI;
    const vector<double>& dt = *DT;
    vector<double> mode(strips);
    vector<map<double, int>> count(strips);

    for (int i = 0; i < bi.size(); i++) {
        count[bi[i]-1][dt[i]]++;        
    }

    for (int i = 0; i < strips; i++) {
        int c = 0;
        double t = 0;
        for (const auto[key, value] : count[i]) {
            if (value > c) {
                c = value;
                t = key;
            }
        }
        mode[i] = t;
    }
    return mode;
}

// takes the mean value of each strip and subtracts it from dt and bt
vector<double> center_dt(vector<double>* DT, vector<double>* BT, vector<int>* BI, int strips) {
    vector<double>& dt = *DT;
    vector<double>& bt = *BT;

    // vector<double> offset = calc_mode(BI, DT, strips);
    vector<double> offset = calc_mean(BI, DT, strips);
    apply_offset(DT, BI, offset);
    apply_offset(BT, BI, offset);
    return offset;
}

// ensures that the file path specified in plot::path exists. also sets the print level for ROOT
void setup() {
    filesystem::create_directories(plot::path);
    gErrorIgnoreLevel = kWarning; // we don't want print outputs
    setup_style(); // use the style defined in ../plot_style.cpp
}

// attempt to fit a gaussian to the input, and, optionally, plot it
template <typename T>
map<string, double> gauss_fit(const vector<T> *input, const vector<double> x_axis = {1000, -500, 500}, bool plot = true, bool save = false) {
    const vector<T> &data = *input;
    TCanvas* canvas = new TCanvas("diff", "diff", 600, 600);
    TH1D *h = new TH1D("1D Histogram", "temporary histogram", int(x_axis[0]), x_axis[1], x_axis[2]); 

    // fill the histogram
    for (int i = 0; i < data.size(); i++) {
        h->Fill(data[i]);
    }

    // perform the actual fit
    TF1* gauss = new TF1("gauss", "gaus", x_axis[1], x_axis[2]);
    h->Fit(gauss, "QR");

    // extract the fit parameters
    // TF1 *fit = h->GetFunction("gaus");
    map<string, double> result = {{"amplitude", gauss->GetParameter(0)}, {"mean", gauss->GetParameter(1)}, {"sigma", gauss->GetParameter(2)}};

    // set up the figure
    if (plot) {
        cout << format("Successfully fit a Gaussian to the data. Mean: %1%, sigma: %2%") % result.at("mean") % result.at("sigma") << endl;
        h->GetXaxis()->SetTitle("Time [ns]");
        h->GetYaxis()->SetTitle("Strip");
        h->GetXaxis()->CenterTitle();
        h->GetYaxis()->CenterTitle();
        h->Draw();
        canvas->SetLogz();
        canvas->SetRightMargin(0.15);
        canvas->Modified(); canvas->Update();
        canvas->WaitPrimitive();
    }
    if (save) { // this is mainly meant for debugging, and currently overwrites the same file if multiple fits are saved. don't use batch mode instead
        string path = plot::path + "gauss_fit" + plot::format;
        canvas->SaveAs(path.c_str());
    }
    canvas->Close();
    h->Delete();
    return result;
}

// merge the analyzed files from output/X.root with the matched files from matched/Xm.root
void merge(int num, char *path[]) {
    for (int i = 1; i < num; i++) { // skip first entry (that's the path to this script)
        // path to analyzed root file
        const char *apath = path[i];

        // path to matched root file
        filesystem::path p(apath);
        string tmp = "match/" + p.stem().string() + "m.root"; //stem is filename without extension; stem(x/y.z) = y
        const char *mpath = tmp.c_str();

        // open the root files
        TFile *fa = TFile::Open(apath);
        TFile *fm = TFile::Open(mpath);
        TTree *ta = (TTree *)fa->Get("tree");
        TTree *tm = (TTree *)fm->Get("a101");

        // define variables from analyzed tree
        Int_t mi[3], mul, N;
        Double_t pt, px[3], py[3], pz[3], deltaE, E_cm[3], E_lab[3], exC12, theta_lab[3], theta_cm[3], phi_lab[3], phi_cm[3];
        ta->SetBranchAddress("N", &N);
        ta->SetBranchAddress("mi", &mi);
        ta->SetBranchAddress("mul", &mul);
        ta->SetBranchAddress("p_tot", &pt);
        ta->SetBranchAddress("px", px);
        ta->SetBranchAddress("py", py);
        ta->SetBranchAddress("pz", pz);
        ta->SetBranchAddress("deltaE", &deltaE);
        ta->SetBranchAddress("E_cm", &E_cm);
        ta->SetBranchAddress("E_lab", &E_lab);
        ta->SetBranchAddress("exC12", &exC12);
        ta->SetBranchAddress("theta_cm", &theta_cm);
        ta->SetBranchAddress("phi_cm", &phi_cm);
        ta->SetBranchAddress("theta_lab", &theta_lab);
        ta->SetBranchAddress("phi_lab", &phi_lab);

        // define variables from matched tree
        // we must be very generous with the allocated space here, since otherwise we risk writing outside the bounds. 100 is probably overkill, though
        UInt_t FT[100], BT[100], FI[100], BI[100], ID[100]; 
        Double_t FE[100], BE[100];
        tm->SetBranchAddress("FI", &FI);
        tm->SetBranchAddress("FT", &FT);
        tm->SetBranchAddress("BI", &BI);
        tm->SetBranchAddress("BT", &BT);
        tm->SetBranchAddress("FE", &FE);
        tm->SetBranchAddress("BE", &BE);
        tm->SetBranchAddress("id", &ID);

        // define destination tree and set its branches
        int ft[3], bt[3], fi[3], bi[3], id[3]; 
        double fe[3], be[3]; // write variables
        string dest = "merged/" + p.filename().string(); // file destination
        TFile f(dest.c_str(), "recreate");
        TTree t("tree", "merged tree for TDC calibration");
        t.Branch("mul", &mul);
        t.Branch("pt", &pt);
        t.Branch("deltaE", &deltaE);
        t.Branch("exC12", &exC12);
        t.Branch("FI", &fi, "FI[3]/I");
        t.Branch("BI", &bi, "BI[3]/I");
        t.Branch("FT", &ft, "FT[3]/I");
        t.Branch("BT", &bt, "BT[3]/I");
        t.Branch("FE", &fe, "FE[3]/D");
        t.Branch("BE", &be, "BE[3]/D");
        t.Branch("id", &id, "id[3]/I");
        t.Branch("px", &px, "px[3]/D");
        t.Branch("py", &py, "py[3]/D");
        t.Branch("pz", &pz, "pz[3]/D");
        t.Branch("E_cm", &E_cm, "E_cm[3]/D");
        t.Branch("E_lab", &E_lab, "E_lab[3]/D");
        t.Branch("theta_cm", &theta_cm, "theta_cm[3]/D");
        t.Branch("phi_cm", &phi_cm, "phi_cm[3]/D");
        t.Branch("theta_lab", &theta_lab, "theta_lab[3]/D");
        t.Branch("phi_lab", &phi_lab, "phi_lab[3]/D");

        // loop over every event in the analyzed tree ta
        for (double i = 0; i < ta->GetEntries(); i++) {
            ta->GetEntry(i); // try to get entry i
            tm->GetEntry(N); // N is the corresponding index in tm for each event in ta
            
            // loop over the three alpha particles
            // we use the multiplicity index mi to get the correct values from the read variables
            for (int i = 0; i < 3; i++) {
                int index = mi[i];
                if (index == -1) { // if it is a reconstructed event
                    ft[i] = 0;
                    bt[i] = 0;
                    fi[i] = -1;
                    bi[i] = -1;
                    fe[i] = 0;
                    be[i] = 0;
                    id[i] = -1;
                } else { // if it is a true event
                    ft[i] = FT[index];
                    bt[i] = BT[index];
                    fi[i] = FI[index];
                    bi[i] = BI[index];
                    fe[i] = FE[index];
                    be[i] = BE[index];
                    id[i] = ID[index];
                }
            }
            t.Fill(); // fill the entries into the merged tree t
        }
        t.Write(); // write to disk
        std::cout << "\033[1;32m" << boost::format("Successfully built %1%") % dest << "\033[0m" << endl; // green colour
    }
}

// write a data_container to disk
void save(data_container* data, string destination) {
    print_title("*** SAVING DATA ***");

    TFile f(destination.c_str(), "recreate"); // open the file
    TTree t("tree", "data from calibrate");

    // define the storage variables for the fill loop
    double pt, deltaE, exC12, px[3], py[3], pz[3], eCM[3], eLab[3], thetaLab[3], phiLab[3], dt[3];
    // double thetaCM[3], phiCM[3];
    int mul, id[3];

    // define the branches in the new tree
    t.Branch("mul", &mul);
    t.Branch("exC12", &exC12);
    t.Branch("p_tot", &pt);
    t.Branch("deltaE", &deltaE); 
    t.Branch("DT", &dt, "dt[3]/D");
    t.Branch("px", &px, "px[3]/D");
    t.Branch("py", &py, "py[3]/D");
    t.Branch("pz", &pz, "pz[3]/D");
    t.Branch("E_cm", &eCM, "E_cm[3]/D");
    t.Branch("E_lab", &eLab, "E_lab[3]/D");
    t.Branch("theta_lab", &thetaLab, "theta_lab[3]/D");
    t.Branch("phi_lab", &phiLab, "phi_lab[3]/D");
    // t.Branch("theta_cm", &thetaCM, "theta_cm[3]/D");
    // t.Branch("phi_cm", &phiCM, "phi_cm[3]/D");
    t.Branch("ID", &id, "ID[3]/i");

    // acquire pointers to the actual data in the data_container
    data->add_dt(); // ensure dt exists
    vector<double> &DT = *data->get_double("dt");
    vector<double> &ECM = *data->get_double("E_cm");
    vector<double> &ELAB = *data->get_double("E_lab");
    vector<double> &PT = *data->get_double("p_tot");
    vector<double> &PX = *data->get_double("px");
    vector<double> &PY = *data->get_double("py");
    vector<double> &PZ = *data->get_double("pz");
    vector<double> &DELTAE = *data->get_double("deltaE");
    vector<double> &THETALAB = *data->get_double("theta_lab");
    vector<double> &PHILAB = *data->get_double("phi_lab");
    // vector<double> &THETACM = *data->get_double("theta_cm");
    // vector<double> &PHICM = *data->get_double("phi_cm");
    vector<double> &EXC12 = *data->get_double("exC12");
    vector<int> &MUL = *data->get_int("mul");
    vector<int> &ID = *data->get_int("ID");

    int n = data->get_n();
    int m = n/3;

    for (int i1 = 0; i1 < m; i1++) {
        int i2 = i1 + m;
        int i3 = i1 + 2*m;

        // single-valued columns
        mul = MUL[i1];
        pt = PT[i1];
        deltaE = DELTAE[i1];
        exC12 = EXC12[i1];

        // vector columns. consider making a simple lambda function which enters all three values to avoid mistakes
        dt[0] = DT[i1]; dt[1] = DT[i2]; dt[2] = DT[i3];
        px[0] = PX[i1]; px[1] = PX[i2]; px[2] = PX[i3];
        py[0] = PY[i1]; py[1] = PY[i2]; py[2] = PY[i3];
        pz[0] = PZ[i1]; pz[1] = PZ[i2]; pz[2] = PZ[i3];
        eCM[0] = ECM[i1]; eCM[1] = ECM[i2]; eCM[2] = ECM[i3];
        eLab[0] = ELAB[i1]; eLab[1] = ELAB[i2]; eLab[2] = ELAB[i3];
        thetaLab[0] = THETALAB[i1]; thetaLab[1] = THETALAB[i2]; thetaLab[2] = THETALAB[i3];
        phiLab[0] = PHILAB[i1]; phiLab[1] = PHILAB[i2]; phiLab[2] = PHILAB[i3];
        // thetaCM[0] = THETACM[i1]; thetaCM[1] = THETACM[i2]; thetaCM[2] = THETACM[i3];
        // phiCM[0] = PHICM[i1]; phiCM[1] = PHICM[i2]; phiCM[2] = PHICM[i3];
        id[0] = ID[i1]; id[1] = ID[i2]; id[2] = ID[i3];

        t.Fill();
    }
    t.Write();
    cout << "Data successfully written to disk as " << destination << "." << endl;
}

// checks if all merged files exists. if even one of them is missing, they are all recreated
// returns the path to each merged file as a string vector
vector<string> check_files(int argc, char *argv[]) {
    vector<string> paths(argc-1); // merged file paths
    bool mergeflag = false;
    std::cout << "\nChecking if all prerequisite files exists..." << endl;
    for (int i = 1; i < argc; i++) {
        filesystem::path p(argv[i]);
        string file = "merged/" + p.filename().string();
        if (!filesystem::exists(file)) {
            mergeflag = true;
            std::cout << "\033[1;31m" << boost::format("Missing %1%") % file << "\033[0m" << endl; // red colour
        }
        paths[i-1] = file; // loop starts at 1 while paths start at 0
    }
    if (mergeflag) {
        std::cout << "Files missing, recreating all merged files..." << endl;
        merge(argc, argv);
    } else {
        std::cout << "\033[1;32m" << "All merged files found." << "\033[0m" << endl; // green colour
    }
    return paths;
}

// flatten the data and extract the necessary branches
// note also the rough filter imposed on it
void prepare_data(int argc, char *argv[], data_container* container, string filter) {
    print_title("*** PREPARING INPUT ***");
    // check if all the merged files exists, and create them if not
    vector<string> paths = check_files(argc, argv);

    std::cout << "Creating a dataframe from the files. This may take a while..." << endl;

    data_container &data = *container;
    vector<double> *FT = new vector<double>(), 
                    *BT = new vector<double>(), 
                    *FE = new vector<double>(), 
                    *BE = new vector<double>(), 
                    *deltaE = new vector<double>(), 
                    *E_cm = new vector<double>(),
                    *E_lab = new vector<double>(),
                    *theta_lab = new vector<double>(),
                    *phi_lab = new vector<double>(),
                    // *theta_cm = new vector<double>(),
                    // *phi_cm = new vector<double>(),
                    *px = new vector<double>(),
                    *py = new vector<double>(),
                    *pz = new vector<double>();
    vector<int> *FI = new vector<int>(), 
                *BI = new vector<int>(), 
                *ID = new vector<int>(); 
    TChain chain("tree");
    for (int i = 0; i < argc-1; i++) {
        chain.Add(paths[i].c_str());
    }
    ROOT::RDataFrame df(chain);
    // auto df_data = df.Filter("mul == 3 && pt < 50e3 && abs(deltaE) < 500"); // energy and momentum conservation
    
    auto df_data = df.Filter(filter); // I think the concentration of data at dt = 0 skews the results here; can maybe be fixed by excluding those
    // auto df_data = df;
    for (int i = 0; i < 3; i++) {
        // extract multi-column data 
        append(FT, df_data.Define("x", (format("FT[%1%]*1e-3") % i).str()).Take<double>("x").GetValue());
        append(BT, df_data.Define("x", (format("BT[%1%]*1e-3") % i).str()).Take<double>("x").GetValue());
        append(FE, df_data.Define("x", (format("FE[%1%]") % i).str()).Take<double>("x").GetValue());
        append(BE, df_data.Define("x", (format("BE[%1%]") % i).str()).Take<double>("x").GetValue());
        append(E_cm, df_data.Define("x", (format("E_cm[%1%]") % i).str()).Take<double>("x").GetValue());
        append(px, df_data.Define("x", (format("pz[%1%]") % i).str()).Take<double>("x").GetValue());
        append(py, df_data.Define("x", (format("py[%1%]") % i).str()).Take<double>("x").GetValue());
        append(pz, df_data.Define("x", (format("pz[%1%]") % i).str()).Take<double>("x").GetValue());
        append(E_lab, df_data.Define("x", (format("E_lab[%1%]") % i).str()).Take<double>("x").GetValue());
        append(theta_lab, df_data.Define("x", (format("theta_lab[%1%]") % i).str()).Take<double>("x").GetValue());
        append(phi_lab, df_data.Define("x", (format("phi_lab[%1%]") % i).str()).Take<double>("x").GetValue());
        // append(theta_cm, df_data.Define("x", (format("theta_cm[%1%]") % i).str()).Take<double>("x").GetValue());
        // append(phi_cm, df_data.Define("x", (format("phi_cm[%1%]") % i).str()).Take<double>("x").GetValue());
        append(FI, df_data.Define("x", (format("FI[%1%]") % i).str()).Take<int>("x").GetValue());
        append(BI, df_data.Define("x", (format("BI[%1%]") % i).str()).Take<int>("x").GetValue());
        append(ID, df_data.Define("x", (format("id[%1%]") % i).str()).Take<int>("x").GetValue());
        // if you need any single-column data, you need to stack it on top of itself thrice which can be done pretty easily in this loop here
        // note that this is only if you need it during this tdc calibration - if you only need it afterwards, simply add it to the dataframe as done below
    }
    
    // add everything to the data container
    data.add_data("FT", *FT);
    data.add_data("BT", *BT);
    data.add_data("FI", *FI);
    data.add_data("BI", *BI);
    data.add_data("ID", *ID);

    // the data below here are not used for anything, and simply tag along for the journey. they will also be filtered by the gauss_cut method
    data.add_data("FE", *FE);
    data.add_data("BE", *BE);
    data.add_data("E_cm", *E_cm);
    data.add_data("px", *px);
    data.add_data("py", *py);
    data.add_data("pz", *pz);
    data.add_data("E_lab", *E_lab);
    data.add_data("theta_lab", *theta_lab);
    data.add_data("phi_lab", *phi_lab);
    // data.add_data("theta_cm", *theta_cm);
    // data.add_data("phi_cm", *phi_cm);
    data.add_data("deltaE", df_data.Take<double>("deltaE").GetValue());
    data.add_data("mul", df_data.Take<int>("mul").GetValue());
    data.add_data("p_tot", df_data.Take<double>("pt").GetValue());
    data.add_data("exC12", df_data.Take<double>("exC12").GetValue());

    std::cout << format("Dataframe was flattened from %1% to %2% entries.") % df_data.Count().GetValue() % FT->size() << endl;
}

void prepare_match(int argc, char *argv[], data_container* container) {
    data_container &data = *container;
    TChain chain("a101");
    for (int i = 1; i < argc; i++) {
        chain.Add(argv[i]);
    }
    ROOT::RDataFrame df(chain);
    auto df_data = df.Filter("mul > 0");
    vector<double> FT = df_data.Define("x", "FT[0]*1e-3").Take<double>("x").GetValue();
    vector<double> BT = df_data.Define("x", "BT[0]*1e-3").Take<double>("x").GetValue();
    vector<double> FE = df_data.Define("x", "FE[0]").Take<double>("x").GetValue();
    vector<double> BE = df_data.Define("x", "BE[0]").Take<double>("x").GetValue();
    auto tmp = df_data.Define("x", "id[0]").Take<uint>("x").GetValue();
    vector<int> ID(tmp.begin(), tmp.end());
    
    // add everything to the data container
    data.add_data("FT", FT);
    data.add_data("BT", BT);
    data.add_data("ID", ID);
    data.add_data("FE", FE);
    data.add_data("BE", BE); 
}