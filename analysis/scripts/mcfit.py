import numpy as np
from matplotlib import pyplot as plt
from iminuit import Minuit
from scipy.linalg import block_diag
from scipy.stats import chi2

class mcfitter:
    def __init__(self,bins,data,transform=None):
        self.bins = bins
        self.data = data
        self.mc_data = None
        self.mc_weight = None
        self.labels = []
        self.transform = lambda x,y: (x,y)
        if transform is not None:
            self.transform = transform

        #Fit results
        self.val = []
        self.cov = []

    def add_mc(self,mc,weight=None,label=None):
        self.labels.append(label)
        if weight is None:
            weight = mc*0+1
        if self.mc_data is None:
            self.mc_data = np.array([mc])
            self.mc_weight = np.array([weight])
        else:
            self.mc_data = np.append(self.mc_data,[mc],axis=0)
            self.mc_weight = np.append(self.mc_weight,[weight],axis=0)

    def prop_error(self,d,cov,df):
        df = np.array([f(d) for f in df])
        return df@cov@df.T

    def fitted_mc(self):
        P,w = self.transform(self.val,self.mc_weight)
        N_D = sum(self.data)
        N = np.sum(w*self.mc_data,axis=1)    # N[j]
        p = N_D*P/N             # p[j] = N_D*P[j]/N_j
        return (p * (w*self.mc_data).T).T
        #return [N_data*self.val[i]*self.mc_weight[i]*self.mc_data[i]/(self.mc_data[i]*self.mc_weight[i]).sum() for i in range(len(self.mc_data))]
        #return [N_data*P*d/d.sum() for P,d in zip(list(self.val),self.mc_data)]

    def prediction(self):
        return sum(self.fitted_mc())
    
    #TODO weight
    def prediction_err(self):
        mc_data = np.array(self.mc_data)
        N_data = sum(self.data)
        N = np.sum(mc_data,axis=1)
        df = [lambda x: [N_data*x[len(N)+i]/N[i] for i in range(len(N))]\
                +[N_data*x[i]/N[i] for i in range(len(N))]]
        fit_err = []

        for i in range(len(mc_data[0])):
            poisson_variance = (np.diag(mc_data.T[i]))
            cov_tot = block_diag(self.cov,poisson_variance)
            data_tot = np.concatenate((self.val,mc_data.T[i]))
            fit_err.append(self.prop_error(data_tot,cov_tot,df)[0][0])
        return np.sqrt(np.array(fit_err))

    def plot_all(self):
        hists = np.array(self.fitted_mc())
        hists = list(np.cumsum(hists,axis=0))
        for i in range(len(hists)-1,-1,-1):
            plt.hist(self.bins[:-1],self.bins,weights=hists[i],\
                    label=self.labels[i]+", "+str(round(self.val[i]*100,2))+"%",linewidth=0.3)
        plt.errorbar(self.bins[:-1]+(self.bins[1]-self.bins[0])/2,self.data,self.data**0.5,\
                color='black',ls='',marker='.',label="Data")
        #plt.errorbar(self.bins[:-1]+(self.bins[1]-self.bins[0])/2,self.prediction(),self.prediction_err(),\
        #        color='red',ls='',marker='',label="Prediction error")
    
    def report_fit(self):
        print("Fit result:")
        for i in range(len(self.val)):
            print("\tp"+str(i)+":",str(round(self.val[i]*100,2))\
                    +"("+str(round(100*np.sqrt(self.cov[0][0]),2))+") %")
        print("\nCorrelation matrix:")
        diag = np.sqrt(np.diag(self.cov))
        cor = self.cov/np.outer(diag,diag)
        s = np.array2string(cor,precision=3).split("\n")
        for s in s: print("\t"+s)
        print("chi2:",self.chisquare(),"/",self.ndf())

    def fit(self,start=None,limit=None,error=None,report=True,no_weight_mode=False):
        if no_weight_mode:
            self.mc_data = self.mc_data*self.mc_weight
            self.mc_weight = self.mc_weight*0+1
        mc_data = self.mc_data
        mc_weight = self.mc_weight
        data = self.data

        def f(p):
            P,w = self.transform(p,mc_weight)
            if sum(P<=0)>0: print("ERROR -> negative p-value:",P)
            return self.lnL2(data,mc_data,P,w)

        l = len(self.mc_data)

        if start is None: start =[1/(l)]*l
        if limit is None: limit = [(0,1)]*l
        if error is None: error = [0.3]*l

        m = Minuit.from_array_func(f,start=start,limit=limit,errordef=1,error=error)
        #m.set_print_level(1)
        m.migrad()
        m.hesse()
        

        self.val = m.np_values()
        self.cov = m.np_covariance()

        if report: self.report_fit()
        return self.val,self.cov
    
    def wrongfit(self,start=None,limit=None,error=None,report=True,no_weight_mode=False):
        if no_weight_mode:
            self.mc_data = self.mc_data*self.mc_weight
            self.mc_weight = self.mc_weight*0+1
        mc_data = self.mc_data
        mc_weight = self.mc_weight
        data = self.data

        def f(p):
            print(p)
            P,w = self.transform(p,mc_weight)
            if sum(P<=0)>0: print("ERROR -> negative p-value:",P)
            
            a = mc_data
            d = data
            N_D = sum(d)
            N = np.sum(w*a,axis=1)    # N[j]
            p = N_D*P/N             # p[j] = N_D*P[j]/N_j

            f = (p@(w*a))

            r = d*np.log(f) - f # sum(a*log(A)-A)
            r[f==0] = 0
            r[np.isnan(r)] = 0
            r = np.sum(r)
            print(r)
            return r
            return self.lnL2(data,mc_data,P,w)

        l = len(self.mc_data)

        if start is None: start =[1/(l)]*l
        if limit is None: limit = [(0,1)]*l
        if error is None: error = [0.3]*l

        m = Minuit.from_array_func(f,start=start,limit=limit,errordef=1,error=error)
        #m.set_print_level(1)
        m.migrad()
        m.hesse()
        

        self.val = m.np_values()
        self.cov = m.np_covariance()

        if report: self.report_fit()
        return self.val,self.cov

    def lnL2(self,d,a,P,w):
        R = np.zeros(d.shape[0])
        index = np.arange(R.shape[0])

        N_D = sum(d)
        N = np.sum(w*a,axis=1)    # N[j]
        p = N_D*P/N             # p[j] = N_D*P[j]/N_j

    # d[i] = 0 case
        # sum(d[i]*log(f[i])-f[i])+ sumsum(a[j,i]*log(A[j,i])-A[j,i])
        # ->
        # -sum(f[i]) + sumsum(a[j,i]*log(A[j,i])-A[j,i])
        d_iszero = (d==0)
        a_ = a[:,d_iszero]
        w_ = w[:,d_iszero]

        t = 1
        A = (a_.T/(1+p*w_.T*t)).T

        f = (p@(w_*A))

        nonzero = a_!=0
        nonzero = np.logical_and(A!=0,a_!=0)
        r = -np.sum(f) + np.sum(a_[nonzero]*np.log(A[nonzero])-A[nonzero]) # sum(a*log(A)-A)

        # BINWISE
        index_ = index[d_iszero]
        temp = a_*np.log(A)-A
        temp[np.isnan(temp)] = 0
        R[index_] = -f + temp.sum(axis=0) # sum(a*log(A)-A)


    # Empty mc case
        wp = (p*w.T).T
        a_temp = a.copy()
        a_temp[wp!=wp.max(axis=0)]=0

        a_case = a_temp.sum(axis=0)==0
        c = np.logical_and(~d_iszero,a_case)
        a_ = a[:,c]
        w_ = w[:,c]
        d_ = d[c]
        A = a_*0.0

        maxpw = np.max((p*w_.T).T,axis=0) # Maximum pw
        pwmax = (p*w_.T).T==maxpw # Where pw is max
        # Calculate A_ki
        a_temp = a_.copy()
        a_temp[pwmax] = 0
        wp_temp = (p*w_.T).T
        wp_temp[pwmax] = 0
        
        k = a_temp*wp_temp/(maxpw-wp_temp)
        k[pwmax] = 0
        k = (d_/(1.0+maxpw) - k.sum(axis=0))/np.sum((p*w_.T).T==maxpw,axis=0)
        
        # TODO do line by line, since pwmax.size != k.size in general
        A[pwmax] = np.tile(k,(A.shape[0],1))[pwmax]
    #print(A)
        
        good = A.sum(axis=0)>0 # Works because A==0 where !pwmax
        a2 = a_[:,~good]
        w2 = w_[:,~good]
        d2 = d_[~good]
        A = A[:,good]
        a_ = a_[:,good]
        w_ = w_[:,good]
        d_ = d_[good]

        #Calculate A_ji
        maxpw = np.max((p*w_.T).T,axis=0) # Maximum pw
        pwmax = (p*w_.T).T==maxpw # Where pw is max
        #print(a_.shape)
        A[~pwmax] = (a_/(1-(p*w_.T).T/((p*w_.T).T).max(axis=0)))[~pwmax]
        #A[~pmax,:] = (a_[~pmax,:]/(1-p[~pmax,np.newaxis]/max(p)))
        f = (p @(w_*A))
    #print("f",f)
        r+=np.sum(d_[f>0]*np.log(f[f>0])-f[f>0])+np.sum(a_[a_>0]*np.log(A[a_>0]))-np.sum(A[A>0])
        
        # BINWISE
        index_ = index[c]
        index2 = index_[~good]
        index_ = index_[good]
        temp = d_*np.log(f)-f
        temp[temp<0] = 0
        R[index_] += temp
        temp = a_*np.log(A)-A
        temp[np.isnan(temp)] = 0
        R[index_] += temp.sum(axis=0)

        #R[index_] = np.sum(d_[f>0]*np.log(f[f>0])-f[f>0])+np.sum(a_[a_>0]*np.log(A[a_>0]))-np.sum(A[A>0])

    # normal case
        special_case = np.sum((a[p==max(p),:]),axis=0)==0
        candidate = np.logical_and(~d_iszero,~a_case)

        a_ = np.concatenate((a[:,candidate],a2),axis=1)
        w_ = np.concatenate((w[:,candidate],w2),axis=1)
        d_ = np.concatenate((d[candidate],d2))
        t = a_[0]*0.0
        
        def err(t): return -d_/(1.0-t) + np.sum((p*(w_*a_).T).T / (1+(p*w_.T).T*t),axis=0)

        positive = -1/np.max((p*w_.T).T,axis=0)
        negative = np.array([1.0]*len(t))

        middle = positive*0.0

        tol = 1e-15 # Tolerance on t, not err(t)!!
        N = int(np.log(tol/(1+1/max(p)))/np.log(0.5))
        for _ in range(N):
            E = err(middle)
            c = E>0
            positive[c] = middle[c]
            negative[~c] = middle[~c]
            middle = (positive+negative)/2
        t = middle
        
        A = a_/(1+(p*w_.T).T*t)
        f = (p @ (w_*A))
        r+=(np.sum(d_*np.log(f)-f)+np.sum(a_[a_!=0]*np.log(A[a_!=0])-A[a_!=0]))
        
        # BINWISE
        index_ = np.concatenate((index[candidate],index2))
        R[index_] = d_*np.log(f)-f
        temp = a_*np.log(A)-A
        temp[np.isnan(temp)] = 0
        R[index_] += temp.sum(axis=0)
        
        r = -2*r
        return r
    
    def chisquare(self,p=None):
        # -2*lnL(y;n)+2*ln(n;n)
        if p is None: p = self.val
        d = self.data
        a = self.mc_data
        mc_weight = self.mc_weight
        P,w = self.transform(p,mc_weight)

        return self.lnL2(d,a,P,w) + 2*(np.sum(d[d>0]*np.log(d[d>0])-d[d>0])+np.sum(a[a>0]*np.log(a[a>0]))-np.sum(a[a>0]))

    def ndf(self):
        return ((self.data+self.mc_data.sum(axis=0))!=0).sum()-len(self.val)

    def pvalue(self,p=None):
        return 1-chi2.cdf(self.chisquare(p=p),self.ndf())

    def mc_pvalue(self,trials,p=None):
        if p is None: p = self.val
        d = self.data
        a = self.mc_data
        mc_weight = self.mc_weight
        P,w = self.transform(p,mc_weight)
        
        lnL = -self.lnL2(d,a,P,w)/2

        N_D = sum(d)
        N = np.sum(w*a,axis=1)    # N[j]
        p = N_D*P/N             # p[j] = N_D*P[j]/N_j
        f = p  @ (w*a)
        
        s = 0.0
        f_ = f.copy()
        a_ = a.copy()
        for _ in range(trials):
            f_ = np.random.multinomial(f.sum(),f/f.sum())
            for i in range(a.shape[0]):
                a_[i,:] = np.random.multinomial(a[i,:].sum(),a[i,:]/a[i,:].sum())
            lnL_ = -self.lnL2(f_,a_,P,w)/2
            if lnL_<=lnL: s+=1
        print(s/trials)
        





    def lnL_binwise(self):
        p = self.val
        d = self.data
        a = self.mc_data
        mc_weight = self.mc_weight
        P,w = self.transform(p,mc_weight)
        # Intro over

        R = np.zeros(d.shape[0])
        index = np.arange(R.shape[0])

        N_D = sum(d)
        N = np.sum(w*a,axis=1)    # N[j]
        p = N_D*P/N             # p[j] = N_D*P[j]/N_j

    # d[i] = 0 case
        # sum(d[i]*log(f[i])-f[i])+ sumsum(a[j,i]*log(A[j,i])-A[j,i])
        # ->
        # -sum(f[i]) + sumsum(a[j,i]*log(A[j,i])-A[j,i])
        d_iszero = (d==0)
        a_ = a[:,d_iszero]
        w_ = w[:,d_iszero]

        t = 1
        A = (a_.T/(1+p*w_.T*t)).T

        f = (p@(w_*A))

        nonzero = a_!=0
        nonzero = np.logical_and(A!=0,a_!=0)
        r = -np.sum(f) + np.sum(a_[nonzero]*np.log(A[nonzero])-A[nonzero]) # sum(a*log(A)-A)

        # BINWISE
        index_ = index[d_iszero]
        temp = a_*np.log(A)-A
        temp[np.isnan(temp)] = 0
        R[index_] = -f + temp.sum(axis=0) # sum(a*log(A)-A)


    # Empty mc case
        wp = (p*w.T).T
        a_temp = a.copy()
        a_temp[wp!=wp.max(axis=0)]=0

        a_case = a_temp.sum(axis=0)==0
        c = np.logical_and(~d_iszero,a_case)
        a_ = a[:,c]
        w_ = w[:,c]
        d_ = d[c]
        A = a_*0.0

        maxpw = np.max((p*w_.T).T,axis=0) # Maximum pw
        pwmax = (p*w_.T).T==maxpw # Where pw is max
        # Calculate A_ki
        a_temp = a_.copy()
        a_temp[pwmax] = 0
        wp_temp = (p*w_.T).T
        wp_temp[pwmax] = 0
        
        k = a_temp*wp_temp/(maxpw-wp_temp)
        k[pwmax] = 0
        k = (d_/(1.0+maxpw) - k.sum(axis=0))/np.sum((p*w_.T).T==maxpw,axis=0)
        
        # TODO do line by line, since pwmax.size != k.size in general
        A[pwmax] = np.tile(k,(A.shape[0],1))[pwmax]
    #print(A)
        
        good = A.sum(axis=0)>0 # Works because A==0 where !pwmax
        a2 = a_[:,~good]
        w2 = w_[:,~good]
        d2 = d_[~good]
        A = A[:,good]
        a_ = a_[:,good]
        w_ = w_[:,good]
        d_ = d_[good]

        #Calculate A_ji
        maxpw = np.max((p*w_.T).T,axis=0) # Maximum pw
        pwmax = (p*w_.T).T==maxpw # Where pw is max
        #print(a_.shape)
        A[~pwmax] = (a_/(1-(p*w_.T).T/((p*w_.T).T).max(axis=0)))[~pwmax]
        #A[~pmax,:] = (a_[~pmax,:]/(1-p[~pmax,np.newaxis]/max(p)))
        f = (p @(w_*A))
    #print("f",f)
        r+=np.sum(d_[f>0]*np.log(f[f>0])-f[f>0])+np.sum(a_[a_>0]*np.log(A[a_>0]))-np.sum(A[A>0])
        
        # BINWISE
        index_ = index[c]
        index2 = index_[~good]
        index_ = index_[good]
        temp = d_*np.log(f)-f
        temp[temp<0] = 0
        R[index_] += temp
        temp = a_*np.log(A)-A
        temp[np.isnan(temp)] = 0
        R[index_] += temp.sum(axis=0)

        #R[index_] = np.sum(d_[f>0]*np.log(f[f>0])-f[f>0])+np.sum(a_[a_>0]*np.log(A[a_>0]))-np.sum(A[A>0])

    # normal case
        special_case = np.sum((a[p==max(p),:]),axis=0)==0
        candidate = np.logical_and(~d_iszero,~a_case)

        a_ = np.concatenate((a[:,candidate],a2),axis=1)
        w_ = np.concatenate((w[:,candidate],w2),axis=1)
        d_ = np.concatenate((d[candidate],d2))
        t = a_[0]*0.0
        
        def err(t): return -d_/(1.0-t) + np.sum((p*(w_*a_).T).T / (1+(p*w_.T).T*t),axis=0)

        positive = -1/np.max((p*w_.T).T,axis=0)
        negative = np.array([1.0]*len(t))

        middle = positive*0.0

        tol = 1e-15 # Tolerance on t, not err(t)!!
        N = int(np.log(tol/(1+1/max(p)))/np.log(0.5))
        for _ in range(N):
            E = err(middle)
            c = E>0
            positive[c] = middle[c]
            negative[~c] = middle[~c]
            middle = (positive+negative)/2
        t = middle
        
        A = a_/(1+(p*w_.T).T*t)
        f = (p @ (w_*A))
        r+=(np.sum(d_*np.log(f)-f)+np.sum(a_[a_!=0]*np.log(A[a_!=0])-A[a_!=0]))
        
        # BINWISE
        index_ = np.concatenate((index[candidate],index2))
        R[index_] = d_*np.log(f)-f
        temp = a_*np.log(A)-A
        temp[np.isnan(temp)] = 0
        R[index_] += temp.sum(axis=0)
        
        return -2*R
    
    def chisquare_binwise(self,p=None):
        # -2*lnL(y;n)+2*ln(n;n)
        if p is None: p = self.val
        d = self.data
        a = self.mc_data
        mc_weight = self.mc_weight
        P,w = self.transform(p,mc_weight)
        dterm = d*np.log(d)
        dterm[d==0] = 0
        aterm = a*np.log(a)
        aterm[a==0] = 0
        return self.lnL_binwise() + 2*(dterm-d+np.sum(aterm,axis=0)-np.sum(a,axis=0))
