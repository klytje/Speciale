

using namespace std;
using namespace ROOT::Math;


class MinuitFitter : public ROOT::Math::IMultiGenFunction {
public:

    MinuitFitter() = default;

    MinuitFitter(vector<Double_t> x, vector<Double_t> y, int ndim, string algoType,
                 std::function<double(double *, const double *)> f) :
            Ndim(ndim), AlgoType(algoType), Xdata(x), Ydata(y), func(f) {
    }


    ROOT::Math::IBaseFunctionMultiDim *Clone() const override {return new MinuitFitter(Xdata, Ydata, Ndim, AlgoType, func);}      
    unsigned int NDim() const override {return Ndim;}
    double getValue(double X, const double *params) const {return func(&X, params);} 


    double DoEval(const double *params) const override {
        Double_t chi = 0;

        // Udregn teoretisk dalitz plot
        for(int xStep = 0; xStep<N; xStep++){
           for(int xStep = 0; xStep<N; xStep++){
              auto x = koordinaten i dit dalitz plot
              auto y = koordinaten i dit dalitz plot
              f(x,y) = ... din function(k, delta, x, y).

              // Udregne chi2 bidraget fra dette bin:
              chi2 += chifunc(fteo, fdata);

           }

        }
        return chi;
    }
}


// Det eneste krav til klassen er at du implementerer:
// Clone, NDim og DoEval. 
// hvor NDim er antallet er fitteparamtre, og DoEval skal returnere din chi2 givet et sæt input parametre
// Du kan så bruge Minimizeren således:
int main(int argc, char const *argv[])
{
std::function<double(double *, const double *)> func;

auto MinFitter{} = new MinuitFitter(Xdata, Ydata, ndim, AlgoType);

MinFitter->setFunction(func);

 auto minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 

auto functor = ROOT::Math::Functor(MinFitter, &MinuitFitter::DoEval, MinFitter->NDim());

minimum->SetFunction(functor);

minimum->Minimize();
}