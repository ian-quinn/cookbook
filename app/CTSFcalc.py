def CalculateCTSF(COND, DEN, CP, Thickness, Resistance, nLTotal):

    COND = [0] + [COND] + [0]
    DEN = [0] + [DEN] + [0]
    CP = [0] + [CP] + [0]
    Thickness = [0] + [Thickness] + [0]
    Resistance = [0.044] + [Resistance] + [0.12]

    X = [0]*125
    XU = [0]*125
    XCV = [0]*125
    LayerAndNodeNumber = [0]*300
    Lamda = [0]*300
    RowCP = [0]*300

    limit = 0.000000001 # limit for response factor calc convergence
    small = 1e-20 # to avoid calculation crash

    # local variables
    nCV = [0]*15


    for nL in range(nLTotal):
        # filt out air resistance layters and assign properties
        if COND[nL] == 0 and not Resistance[nL] == 0:
            COND[nL] = 0.0263
            DEN[nL] = 1.1614
            CP[nL] = 1007
            Thickness[nL] = COND[nL] * Resistance[nL]

        # set 6 as the default nodes of each layer
        nCV[nL] = 6
        # volume thermal properties

    XU[2] = 0
    n2 = 2
    for nL in range(nLTotal):
        nLast = n2
        n1 = nLast + 1
        n2 = nLast + nCV[nL]
        for n in range(n1, n2+1):
            spacing = (n-nLast) / nCV[nL]
            XU[n] = XU[nLast] + Thickness[nL] * spacing
            LayerAndNodeNumber[n-1] = nL

    nMax = n2
    #  of nodes in the constructon
    #  Assign total number of nodes
    nM2 = nMax - 1
    nM3 = nM2 - 1
    X[1] = XU[2]
    for n in range(2, nM2 + 1):
        X[n] = 0.5 *  ( XU[n+1] + XU[n] )
        XCV[n] = XU[n+1] - XU[n]
    X[nMax] = XU[nMax]

    # X, XU, XCV # official returns

    # ReDim something
    a = [0] * (nMax+1)
    d = [0] * (nMax+1)
    B = [0] * (nMax+1)
    C = [0] * (nMax+1)
    phi = [0] * (nMax+1)
    XYPHi = [0] * (nMax+1)
    ZPhi = [0] * (nMax+1)

    for I in range(nMax+1):
        phi[I] = 0

    #   assign layer number to each control volume in the construction
    for inode in range(2, nM2+1):
        numLayer = LayerAndNodeNumber[inode] # check which layer it belongs to
        Lamda[inode] = COND[numLayer]
        RowCP[inode] = CP[numLayer] * DEN[numLayer]

    #   start solving for temperature and response factors for each time step
    XFlux = float('inf')
    YFlux = float('inf')
    ZFlux = float('inf')
    XResponse = [0]*300
    YResponse = [0]*300
    ZResponse = [0]*300

    HoursCount = 1
    #
    TimeCount = 0
    TimeStep = 60
    N_TimeSteps = 9600
    for StepCount in range(1, N_TimeSteps+1):
        if TimeCount >= 1*24*3600:
            TimeStep = 60
        TimeCount = TimeCount + TimeStep
        #   solve the outside and cross response factors
        if abs(XFlux) < limit and abs(YFlux) < limit and TimeCount > 7200:
            pass # 好像是个求解的终止条件
        else:
            # Call(GenerateCoefficients)
            # def ApplyBoundaryConditions():
            # 界定了第一个节点的温度值，基本上说明数组是从1开始的？
            if TimeCount <= 3600:
                phi[1] = TimeCount / 3600
            elif TimeCount > 3600 and TimeCount <= 7200:
                phi[1] = 2 - TimeCount / 3600
            else:
                phi[1] = 0
            phi[nMax] = 0

            # def TDMACoefficientsSetup():
            BETA = 4/3

            for I in range(2, nM2+1):
                B[I] = 0
                a[I] = 0
                d[I] = 0
                C[I] = 0
            #   constant volumetric terms
            for I in range(2, nM2+1):
                APT = RowCP[I] / TimeStep
                C[I] = ( C[I] + APT * phi[I] )  * XCV[I]
                d[I] = APT * XCV[I]
            #   Interior nodes
            for I in range(2, nM3+1):
                a[I] = 2 * Lamda[I] * Lamda[I+1] /  ( XCV[I] * Lamda[I+1] + XCV[I+1] * Lamda[I] ) + small
                B[I+1] = a[I]
            #   left handside boundary condition
            B[2] = BETA *  ( Lamda[2] /  ( 0.5 * XCV[2] ) )  + small
            a[1] = B[2]
            B[1] = ( BETA - 1 ) * a[2]
            a[2] = a[2] + B[1]
            C[2] = C[2] + B[2] * phi[1]
            d[2] = d[2] + B[2]
            B[2] = 1
            #   right hand side boundary condition
            a[nM2] = BETA * ( Lamda[nM2] /  ( 0.5 * XCV[nM2] ) ) + small
            B[nMax] = a[nM2]
            a[nMax] = ( BETA - 1 )  * B[nM2]
            B[nM2] = B[nM2] + a[nMax]
            C[nM2] = C[nM2] + a[nM2] * phi[nMax]
            d[nM2] = d[nM2] + a[nM2]
            a[nM2] = 0
            for I in range(2, nM2+1):
                d[I] = d[I] + a[I] + B[I]

            # def TDMASolver():
            # Decomposition and forward substitution.
            PTX = [0]*300
            QTX = [0]*300

            PTX[1] = 0
            QTX[1] = phi[1]
            for I in range(2, nM2+1):
                Denom = d[I] - PTX[I-1] * B[I]
                PTX[I] = a[I] / Denom
                QTX[I] = ( C[I] + B[I] * QTX[I-1] )  / Denom
            #
            # Backsubstitution.
            for J in range(nM2, 1, -1):
                phi[J] = QTX[J] + PTX[J] * phi[J+1]
            #
            #    compute fluxes at the inside and outside faces due to excitation
            #    at the exterior surface
            XFlux = a[1] *  ( phi[1] - phi[2] )  + B[1] *  ( phi[3] - phi[2] )
            YFlux = - ( B[nMax] *  ( phi[nMax] - phi[nM2] )  + a[nMax] *  ( phi[nM3] - phi[nM2] ) )
            #    compute fluxes at the inside surface due to excitation
            #    at the interior surface
            ZFlux = a[1] *  ( phi[1] - phi[2] )  + B[1] *  ( phi[3] - phi[2] )

        #    save the temperatures for external exitation boundary condition
        if TimeCount == HoursCount * 3600:
            XResponse[HoursCount] = XFlux
            YResponse[HoursCount] = YFlux
            HoursCount = HoursCount + 1

    PF = 24
    #  calculate periodic response factor for each construction in the zone
    #  from the surface reposne factors

    XPRF = [0] * 24
    YPRF = [0] * 24
    ZPRF = [0] * 24
    for I in range(24):
        XPRF[I] = XResponse[I+1]
        YPRF[I] = YResponse[I+1]
        ZPRF[I] = ZResponse[I+1]
        for J in range(1, 8):
            if not abs(XResponse[I + J * 24]) == 0:
                XPRF[I] = XPRF[I] + XResponse[I + J * 24]
            if not abs(YResponse[I + J * 24]) == 0:
                YPRF[I] = YPRF[I] + YResponse[I + J * 24]
            if not abs(ZResponse[I + J * 24]) == 0:
                ZPRF[I] = ZPRF[I] + ZResponse[I + J * 24]

    def sumArray(array):
        sigma = 0
        for value in array:
            sigma = sigma + value
        return sigma

    CTSOut = []
    for J in range(24):
        CTSOut.append(YPRF[J] / sumArray(YPRF))

    return CTSOut
