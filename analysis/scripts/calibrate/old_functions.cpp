// optimization stuff
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

// other root stuff
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TLine.h>
#include <TChain.h>
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TROOT.h>

// other stuff
#include <filesystem>
#include <boost/format.hpp>
#include <iostream>
#include <cmath>

// my own stuff
#include "plots.cpp"
#include "utility.cpp"
#include "data_fix.cpp"
#include "calibrate.cpp" // breaks any attempt at compilation. this is not supposed to be compiled anyway

using namespace std;
using boost::format;

void apply_offset(data_container* data, const vector<vector<double>> offset_left, const vector<vector<double>> offset_right, double center) {
    data->add_dt();
    const vector<double> &dt = *data->get_double("dt");
    const vector<int> &fi = *data->get_int("FI");
    const vector<int> &bi = *data->get_int("BI");
    const vector<int> &id = *data->get_int("ID");
    vector<double> &bt = *data->get_double("BT");
    
    for (int i = 0; i < dt.size(); i++) {
        if (center < dt[i]) {
            bt[i] -= offset_right[id[i]][bi[i]-1];
        } else {
            bt[i] -= offset_left[id[i]][bi[i]-1];
        }
    }
    return;
}

// subtracts the offsets specified in offset_left and offset_right from BT depending on whether an element is on the left or right of center
// this method applies the offsets on data from all detectors simultaneously
void apply_offset_old(const vector<double> *FT, vector<double> *BT, const vector<int> *FI, const vector<int> *BI, const vector<int> *ID, const vector<vector<double>> offset_left, const vector<vector<double>> offset_right, double center) {
    const vector<double> &FTref = *FT;
    const vector<int> &FIref = *FI;
    const vector<int> &BIref = *BI;
    const vector<int> &IDref = *ID;
    vector<double> &BTref = *BT;

    vector<double> dt = calc_dt(FT, BT);
    int n = BTref.size();

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

        vector<int> bi;
        // extract the other necessary stuff from the coincidences, and calculate dt
        tie(ft, bt, std::ignore, bi, std::ignore) = find_coincidences(FT, BT, FI, BI, ID, det);
        dt_S3 = calc_dt(&ft, &bt);

        hist_info info(plot::x_axes::standard, plot::y_axis.at(det), "tmp", "Title", "dt [ns]", "Back strip id");
        hist2D(&dt_S3, &bi, info);

        // correct each dt entry in the main dt vector
        for (int i = 0; i < index.size(); i++) {
            cout << "dt before: " << dt[index[i]] << ", dt after: " << dt_S3[i] << endl;
            dt[index[i]] = dt_S3[i];
        }
    }

    //*** W1 detectors ***//
    // for these we apply the offset to any event observed by them
    int id;
    for (int i = 0; i < BT->size(); i++) {
        id = IDref[i];
            if (center < dt[i]) {
                BTref[i] -= offset_right[id][BIref[i]-1];
            } else {
                BTref[i] -= offset_left[id][BIref[i]-1];
            }
    }

    //*** S3 detectors ***//
    // for these we should only apply the offset to the actual events involved in determining it in the first place
    for (int det = 0; det <= 1; det++) { // iterate over the two S3 detectors
        vector<int> index(BT->size());

        // set up the index vector
        for (int i = 0; i < index.size(); i++) {
            index[i] = i;
        }

        // extract the back indices of each coincidence event
        tie(std::ignore, index) = find_coincidences(&index, FI, BI, ID, det);

        vector<double> ft, bt, dt;
        vector<int> bi;
        // extract the other necessary stuff from the coincidences, and calculate dt
        tie(ft, bt, std::ignore, bi, std::ignore) = find_coincidences(FT, BT, FI, BI, ID, det);
        dt = calc_dt(&ft, &bt);
        hist_info info(plot::x_axes::standard, plot::y_axis.at(det), "tmp", "Title", "dt [ns]", "Back strip id");
        hist2D(&dt, &bi, info);

        for (int i = 0; i < index.size(); i++) {
            if (center < bt[i] - ft[i]) bt[i] -= offset_right[det][bi[i]-1];
            else bt[i] -= offset_left[det][bi[i]-1];
            BTref[index[i]] = bt[i];
        }

        tie(ft, bt, std::ignore, bi, std::ignore) = find_coincidences(FT, BT, FI, BI, ID, det);
        dt = calc_dt(&ft, &bt);

        hist2D(&dt, &bi, info);
    }
    return;
}

// a double filter to save an array iteration
template<typename T>
vector<T> detector_dt_filter(vector<T> *input, vector<double> *DT, vector<int> *ID, int detector, double a, double b) {
    int n = ID->size(); 
    vector<int> &IDref = *ID; // this construction does *not* copy the vector, only dereferences it
    vector<double> &DTref = *DT;
    vector<T> &inputref = *input;
    vector<T> result(n);
    int m = 0; // current index in result
    double dt;
    for (int i = 0; i < n; i++) {
        if (IDref[i] == detector) {
            dt = DTref[i];
            if (a < dt && dt < b) {
                result[m] = inputref[i];
                m++;
            }
        }
    }
    result.resize(m);
    return result; // remove all the unset entries    
}

// splits the input into two halves on each side of center, bounded by a and b
// note that it imposes the same dt_filter as the method above
template <typename T>
pair<vector<T>, vector<T>> dt_split(const vector<T> *input, const vector<double> *DT, double center, double a, double b) {
    int n = DT->size();
    const vector<double> &DTref = *DT;
    const vector<T> inputref = *input;
    vector<T> result_left(n);
    vector<T> result_right(n);
    int m_left = 0; 
    int m_right = 0; 
    double dt;
    for (int i = 0; i < n; i++) {
        dt = DTref[i];
        if (center < dt && dt < b) { // right half, center < dt < b
            result_right[m_right] = inputref[i];
            m_right++;
        }
        else if (dt <= center && a < dt) { // left half, a < dt <= center. note <= so we do not lose any values
            result_left[m_left] = inputref[i];
            m_left++;
        }
    }
    result_left.resize(m_left);
    result_right.resize(m_right);
    return pair(result_left, result_right); // remove all the unset entries
}

// since my data is weird and shows *two* columns, one on either side of 0, this method attempts to center them both around 0
// note also how it overwrites dt and bi. it is encapsulated such that everything but dt and bi are cleared afterwards
// I see no reason why it should not work with normal, single-column data either - then one side would just be empty
pair<vector<double>, vector<double>> center_columns(vector<double> *FT, vector<double> *BT, vector<int> *FI, vector<int> *BI, vector<double> *DT, vector<int> *ID, vector<double> x_axis, int strips) {
    int dt_a = x_axis[1], dt_b = x_axis[2]; // the bounds of my two columns
    // split the data into the two columns seen on plot_raw
    auto[bi_left, bi_right] = dt_split(BI, DT, 0, dt_a, dt_b);
    auto[dt_left, dt_right] = dt_split(DT, DT, 0, dt_a, dt_b);

    // calculate the offset for each column
    vector<double> offset_left = calc_offset(&bi_left, &dt_left, strips);
    vector<double> offset_right = calc_offset(&bi_right, &dt_right, strips);

    // apply the offset to the data
    apply_offset(&dt_left, &bi_left, offset_left);
    apply_offset(&dt_right, &bi_right, offset_right);

    // merge the two halves into one
    append(&dt_left, &dt_right);
    // append(&bt_left, &bt_right);

    // before we assign these two halves back into dt, we need to sort the other vectors (otherwise their lengths will differ)
    tie(*FT, *BT, *FI, *BI, *ID) = dt_filter(FT, BT, DT, FI, BI, ID, dt_a, dt_b);

    // the following applies the offset to "bt", such that if we calculate dt = (bt - ft) later, it will correctly be centered
    apply_offset(FT, BT, BI, offset_left, offset_right, 0);
    *DT = dt_left;
    return pair(offset_left, offset_right);
}

// subtracts the offsets specified in offset_left and offset_right from BT depending on whether an element is on the left or right of center
// this correctly centers dt = (bt - ft) when it is calculated afterwards
void apply_offset(const vector<double> *FT, vector<double> *BT, const vector<int> *BI, const vector<double> offset_left, const vector<double> offset_right, double center) {
    const vector<int> &BIref = *BI;
    const vector<double> &FTref = *FT;
    vector<double> &BTref = *BT;
    for (int i = 0; i < BT->size(); i++) {
        if (center < BTref[i] - FTref[i]) {
            BTref[i] -= offset_right[BIref[i]-1];
        } else {
            BTref[i] -= offset_left[BIref[i]-1];
        }
    }
    return;
}

//*** FROM OLD DATA_FIX ***//
template <typename T>
vector<T> get_highest(vector<T>* t, vector<double>* cmp) {
    vector<T> &tref = *t;
    vector<double> &cmpref = *cmp;

    int i2, i3, n = cmpref.size();
    int n_red = n/3;
    double m;
    vector<T> result(n_red);
    for (int i1 = 0; i1 < n_red; i1++) {
        i2 = i1 + n_red; i3 = i2 + n_red;
        m = std::max({cmpref[i1], cmpref[i2], cmpref[i3]});
        if (m == cmpref[i1]) result[i1] = tref[i1];
        else if (m == cmpref[i2]) result[i1] = tref[i2];
        else if (m == cmpref[i3]) result[i1] = tref[i3];
    }
    return result;
}

template <typename T>
vector<T> get_lowest(vector<T>* t, vector<double>* cmp) {
    vector<T> &tref = *t;
    vector<double> &cmpref = *cmp;

    int i2, i3, n = cmpref.size();
    int n_red = n/3;
    double m;
    vector<T> result(n_red);
    for (int i1 = 0; i1 < n_red; i1++) {
        i2 = i1 + n_red; i3 = i2 + n_red;
        m = std::min({cmpref[i1], cmpref[i2], cmpref[i3]});
        if (m == cmpref[i1]) result[i1] = tref[i1];
        else if (m == cmpref[i2]) result[i1] = tref[i2];
        else if (m == cmpref[i3]) result[i1] = tref[i3];
    }
    return result;
}

template <typename T>
vector<T> get_between(vector<T>* t, vector<double>* cmp) {
    vector<T> &tref = *t;
    vector<double> &cmpref = *cmp;

    int i2, i3, n = cmpref.size();
    int n_red = n/3;
    double ma, mi;
    vector<T> result(n_red);
    for (int i1 = 0; i1 < n_red; i1++) {
        i2 = i1 + n_red; i3 = i2 + n_red;
        mi = std::min({cmpref[i1], cmpref[i2], cmpref[i3]});
        ma = std::max({cmpref[i1], cmpref[i2], cmpref[i3]});
        if (mi == cmpref[i2] && ma == cmpref[i3]) result[i1] = tref[i1];
        else if (mi == cmpref[i3] && ma == cmpref[i2]) result[i1] = tref[i1];

        else if (mi == cmpref[i1] && ma == cmpref[i2]) result[i1] = tref[i3];
        else if (mi == cmpref[i2] && ma == cmpref[i1]) result[i1] = tref[i3];

        else if (mi == cmpref[i1] && ma == cmpref[i3]) result[i1] = tref[i2];
        else if (mi == cmpref[i3] && ma == cmpref[i1]) result[i1] = tref[i2];
    }
    return result;
}

// filter the input to elements satisfying a < dE < b
template <typename T>
vector<T> dE_filter(const vector<T> *input, const vector<double> *dE, double a, double b) {
    int n = dE->size();
    const vector<double> &dEref = *dE;
    const vector<T> &inputref = *input;
    vector<T> result(3*n);
    vector<bool> filter(3*n, false);

    int i2, i3;
    for (int i1 = 0; i1 < n; i1++) {
        i2 = i1 + n; i3 = i1 + 2*n;
        if (a < dEref[i1] && dEref[i1] < b) {
            filter[i1] = true;
            filter[i2] = true;
            filter[i3] = true;
        }
    }

    int m = 0; // current index in result
    for (int i = 0; i < 3*n; i++) {
        if (filter[i]) {
            result[m] = inputref[i];
            m++;
        }
    }

    result.resize(m); // remove all the unset entries
    return result;
}

template<typename T>
void check_col(vector<double>* t, vector<T>* cmp, vector<int> *id) {
    vector<double> x = {1000, 65000, 66000};
    hist_info info(x, "title", "tmp", "t [ns]", "count");
    auto id_min = get_lowest(id, cmp);
    auto _min = get_lowest(t, cmp);

    auto id_max = get_highest(id, cmp);
    auto _max = get_highest(t, cmp);

    auto id_between = get_between(id, cmp);
    auto _between = get_between(t, cmp);

    // cout << "plotting min entries" << endl;
    // hist1D(&_min, info);
    // cout << "plotting middle entries" << endl;
    // hist1D(&_between, info);
    // cout << "plotting max entries" << endl;
    // hist1D(&_max, info);

    if (true) { // all data from each detector
    for (int det = 0; det < 4; det++) {
        std::cout << "\033[1;34m" << "\nPlotting for detector " << detector_map.at(det) << "." << "\033[0m" << endl; 
        auto t_min = detector_filter(&_min, &id_min, det);
        auto t_max = detector_filter(&_max, &id_max, det);
        auto t_between = detector_filter(&_between, &id_between, det);
        std::cout << "Plotting max entries" << endl;
        hist1D(&t_max, info);

        std::cout << "Plotting middle entries" << endl;
        hist1D(&t_between, info);

        std::cout << "Plotting min entries" << endl;
        hist1D(&t_min, info);
    }
    }

    if (false) { // multi-detector data only
    std::cout << "\nFiltering to only multi-detector data and repeating the plots." << endl;
    vector<int> &idref = *id;
    vector<double> &tref = *t;
    vector<T> &cmpref = *cmp;
    for (int det = 0; det < 4; det++) {
        std::cout << "\033[1;34m" << "\nPlotting for detector " << detector_map.at(det) << "." << "\033[0m" << endl; 

        int n = t->size();
        int m = n/3;

        vector<double> _t(n);
        vector<T> _cmp(n);
        vector<int> _id(n);
        vector<bool> filter(n, false);

        int i2, i3;
        for (int i1 = 0; i1 < m; i1++) {
            i2 = i1 + m; i3 = i2 + m;
            if (!(idref[i1] == idref[i2] && idref[i1] == idref[i3])) {
                filter[i1] = true;
                filter[i2] = true;
                filter[i3] = true;
            }
        }

        int c = 0; // counter
        for (int i = 0; i < n; i++) {
            if (filter[i]) {
                _t[c] = tref[i];
                _cmp[c] = cmpref[i];
                _id[c] = idref[i];
                c++;
            }
        }
        _t.resize(c); _cmp.resize(c); _id.resize(c);

        id_min = get_lowest(&_id, &_cmp);
        _min = get_lowest(&_t, &_cmp);

        id_max = get_highest(&_id, &_cmp);
        _max = get_highest(&_t, &_cmp);

        id_between = get_between(&_id, &_cmp);
        _between = get_between(&_t, &_cmp);

        auto t_min = detector_filter(&_min, &id_min, det);
        auto t_max = detector_filter(&_max, &id_max, det);
        auto t_between = detector_filter(&_between, &id_between, det);

        std::cout << "Plotting max entries" << endl;
        hist1D(&t_max, info);

        std::cout << "Plotting middle entries" << endl;
        hist1D(&t_between, info);

        std::cout << "Plotting min entries" << endl;
        hist1D(&t_min, info);
    }
    }
}

void debug(data_container* data) {
    vector<double>* FT = data->get_double("FT");
    vector<double>* BT = data->get_double("BT");
    vector<double>* FE = data->get_double("FE");
    vector<double>* BE = data->get_double("BE");
    vector<double>* dE = data->get_double("deltaE");
    vector<int>* ID = data->get_int("ID");
    vector<double> ft, bt, fe, be, id;

    // print_title("Plotting BT filtered on BT");
    // check_col(BT, BT, ID);

    print_title("Plotting BT filtered on FE");
    check_col(BT, FE, ID);

    print_title("Plotting FT filtered on FE");
    check_col(FT, FE, ID);

    cout << "Plotting FT filtered on deltaE" << endl;
    vector<double> x = {1000, 65000, 66000};
    hist_info info(x, "title", "tmp", "t [ns]", "count");

    cout << "deltaE range: -20, 20" << endl;
    ft = dE_filter(FT, dE, -20, 20);
    hist1D(&ft, info);
    cout << "deltaE range: -40, 40" << endl;
    ft = dE_filter(FT, dE, -40, 40);
    hist1D(&ft, info);
    cout << "deltaE range: -60, 60" << endl;
    ft = dE_filter(FT, dE, -60, 60);
    hist1D(&ft, info);
    cout << "deltaE range: -80, 80" << endl;
    ft = dE_filter(FT, dE, -80, 80);
    hist1D(&ft, info);
    cout << "deltaE range: -100, 100" << endl;
    ft = dE_filter(FT, dE, -100, 100);
    hist1D(&ft, info); // missing a lot of entries compared to the next one    
}
