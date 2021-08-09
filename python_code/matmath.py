from scipy import interpolate, optimize


class MatMath:
    @staticmethod
    def interp1(y_arr, chord_arr, y, kind):
        # theta = Mat.interp1(y_arr, theta_arr, y, 'linear')
        fnct = interpolate.interp1d(y_arr, chord_arr, kind=kind)
        chord = fnct(y)
        return chord

    @staticmethod
    def interp1_fnct(y_arr, chord_arr, kind):
        # theta = Mat.interp1(y_arr, theta_arr, y, 'linear')
        fnct = interpolate.interp1d(y_arr, chord_arr, kind=kind)
        return fnct

    @staticmethod
    def interp2(X, Y, V, Xq, Yq, kind):
        # Cl = Mat.interp2(X, Y, V, Xq, Yq, 'linear')
        pass

    @staticmethod
    def fsolve(funzero, xguess, *args):
        # options = optimoptions(
        #           @fsolve, 'Display', 'off', 'MaxFunctionEvaluations', 1500)
        # xguess = lambda_c + 0.1
        # [x0, fval, exitflag, output] = fsolve( @funzero, xguess, options)

        # https://docs.scipy.org/
        # doc/scipy/reference/generated/scipy.optimize.fsolve.html

        # scipy.optimize.fsolve(
        #     func, x0, args=(), fprime=None, full_output=0, col_deriv=0,
        #     xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100,
        #     diag=None)

        x, infodict, ier, mesg = optimize.fsolve(
            funzero, xguess, args=args, fprime=None, full_output=True,
            xtol=1.0e-05)

        x0 = x
        fval = infodict['fvec']
        exitflag = ier
        output = mesg
        return [x0, fval, exitflag, output]
