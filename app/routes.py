# import os, re, json
import math
import datetime
import os

from flask import flash, render_template, jsonify, request
# from werkzeug.urls import url_parse
# from werkzeug.utils import secure_filename

import numpy as np
import pandas as pd

from pyecharts import options as opts
from pyecharts.charts import Bar, Line

# ------------from current app-------------
from app import app
from app.forms import OffsetCalcForm
from app.sunpathcalc import getSEA, getAZ
from app.CTSFcalc import CalculateCTSF

def offsetCalc(city, date, psi, cloudiness, maxTemp, minTemp, Uwall, Uglazing, ampPeople, ampLight, ampEquip, setTemp, SHGC, ctsf, rts, H, Hsill, WWR, cond, den, cp, thickness):
    def loopSigma(inputseries, timeseries):
        def loopBackwards(inputseries, breakpoint):
            if breakpoint > len(inputseries) or breakpoint < 0:
                print('breakpoint out of the series index')
                return inputseries
            return inputseries[breakpoint::-1] + inputseries[-1:breakpoint:-1]

        def dotProduct(list1, list2):
            product = 0
            for i in range(len(list1)):
                product += list1[i]*list2[i]
            return product

        inputseries = inputseries.tolist()
        timeseries = timeseries.tolist()
        if len(inputseries) != len(timeseries):
            print('data and series are not matched')
            return [0 for i in range(len(timeseries))]
        else:
            sigma = []
            for i in range(len(inputseries)):
                recseries = loopBackwards(inputseries, i)
                sigma.append(dotProduct(recseries, timeseries))
            return np.array(sigma)
    # Calculate distance
    def arraySigma(array):
        sigma = 0
        for i in range(len(array)):
            sigma += array[i]
        return sigma

    def heavisideFilter(array):
        if isinstance(array, list):
            for i in range(len(array)):
                if array[i] < 0:
                    array[i] = 0
            return array
        else:
            if array > 0:
                return array
            else:
                return 0

        

    ###########
    # PRESETS #
    ###########
    schPeople = np.array([0,0,0,0,0,0,0,0,0.2,0.8,1,1,0.5,0.5,1,1,1,0.5,0.2,0,0,0,0,0])
    schLight = np.array([0,0,0,0,0,0,0,0,0.5,1,1,1,0.75,0.75,1,1,1,1,0.5,0,0,0,0,0])
    schEquip = np.array([0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.5,0.8,1,1,0.9,0.9,1,1,1,1,0.8,0.5,0.1,0.1,0.1,0.1])
    CTSFdict = {
        'Blink':(0.99,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), \
        'Curtain':(0.18,0.58,0.2,0.04,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), \
        'Stud':(0.06,0.42,0.33,0.13,0.04,0.01,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), \
        'EIFS':(0.02,0.25,0.31,0.2,0.11,0.05,0.03,0.02,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), \
        'Brick':(0,0.05,0.14,0.17,0.15,0.12,0.09,0.07,0.05,0.04,0.03,0.02,0.02,0.01,0.01,0.01,0.01,0.01,00,0,0,0,0,0), \
        'BrickMass':(0.03,0.03,0.03,0.04,0.04,0.04,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.03,0.03)
    }

    RTSsolardict = {
        'Blink':(0.99,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), \
        'Light':(0.45,0.2,0.11,0.07,0.05,0.03,0.02,0.02,0.01,0.01,0.01,0.01,0.01,0,0,0,0,0,0,0,0,0,0,0), \
        'Medium':(0.29,0.15,0.1,0.07,0.06,0.05,0.03,0.03,0.03,0.03,0.02,0.02,0.02,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0,0), \
        'Heavy':(0.27,0.13,0.07,0.05,0.04,0.04,0.03,0.03,0.03,0.03,0.03,0.03,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.01,0.01)
    }

    RTSmassdict = {
        'Blink':(0.99,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), \
        'Light':(0.43,0.19,0.11,0.07,0.05,0.03,0.03,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0,0,0,0,0,0,0,0,0), \
        'Medium':(0.33,0.16,0.1,0.07,0.05,0.04,0.03,0.03,0.02,0.02,0.02,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0), \
        'Heavy':(0.25,0.09,0.06,0.05,0.05,0.04,0.04,0.04,0.03,0.03,0.03,0.03,0.03,0.03,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02)
    }

    # Params formation
    citydict = {'Shanghai':(121.45,31.4), 'Harbin':(126.77,45.75), 'Beijing':(116.47,39.80), 'Guangzhou':(113.33,23.17), 'Lhasa':(91.13,29.67), 'Urumqi':(87.65,43.78), 'Lanzhou':(103.88,36.5), 'Kunming':(102.68,25.02)}
    longitude = citydict[city][0]
    latitude = citydict[city][1]

    timezone = 8

    date = date# if date else '2020-12-21'
    psi = psi# if psi else 180
    # outside srf resistance 0.044 / inside verticle srf resistance 0.120
    Uwall = 1/(1/Uwall+0.164)
    # ...

    RTSsolar = np.array(RTSsolardict[rts])
    RTSmass = np.array(RTSmassdict[rts])
    if ctsf == 'Customize!':
        CTSFwall = np.array(CalculateCTSF(cond, den, cp, thickness, cond/thickness, 3))
    else:
        CTSFwall = np.array(CTSFdict[ctsf])

    # import epw information
    epwpath = os.path.join(app.config['EPW_PATH'], city + '.epw')
    checkdate = datetime.datetime.strptime(date, '%Y-%m-%d')
    month = int(datetime.datetime.strftime(checkdate, '%m'))
    day = int(datetime.datetime.strftime(checkdate, '%d'))
    num_day = int(datetime.datetime.strftime(checkdate, '%j'))
    
    # whether to use sinusoid temperature?
    if maxTemp and minTemp:
        seriesTout = np.array([(maxTemp-minTemp)*math.sin((i-7)*math.pi/12)/2+(maxTemp+minTemp)/2 for i in range(24)])
    else:
        Tout = []
        epw = open(epwpath.strip('\n'), 'r')
        for line in epw:
            checklist = line.split(',')
            if checklist[1] == str(month) and checklist[2] == str(day) and len(checklist) == 35:
                # len(checklist) is to make sure that you locate the right line not the file header
                Tout.append(float(checklist[6]))
        seriesTout = np.array(Tout)
        epw.close()
    
    # Calculate load from solar beam
    azimuths = [getAZ(latitude, longitude, timezone, date+' '+str(i)+':00') for i in range(24)]
    altitudes = [getSEA(latitude, longitude, timezone, date+' '+str(i)+':00') for i in range(24)]
    theta_cos = []
    for i in range(24):
        # if altitudes[i] > 0:
        beta_rad = math.radians(altitudes[i])
        phi_rad = math.radians(abs(azimuths[i]-psi)) #if abs(azimuths[i]-psi) < 90 else math.pi/2
        theta_cos.append(math.cos(beta_rad)*math.cos(phi_rad))
        # else:
            # theta_cos.append(0)


    # whether to use epw recorded radiation?
    Ebeam_raw, Ediffuse_raw = [], []
    if cloudiness or cloudiness == 0:
        Eterrestrial, ab, ad = 1367, 0.64967, 0.2044
        tb = 0.38
        td = 2.5
        # ab = 1.454 - 0.406*tb - 0.268*td + 0.021*tb*td
        # ad = 0.507 + 0.205*tb - 0.08*td - 0.19 *tb*td
        E0 = Eterrestrial * (1 + 0.033 * math.cos((num_day-3)*360/365))
        for i in range(len(altitudes)):
            islit = np.heaviside(altitudes[i],1)
            m = 1 / ( math.sin(math.radians(islit * altitudes[i])) + 0.50572*math.pow((6.07995+math.radians(islit * altitudes[i])),-1.6364) )
            Ebeam_raw.append(E0 * math.exp(-tb*math.pow(m,ab)) * (1-cloudiness) * islit)
            Ediffuse_raw.append(E0 * math.exp(-td*math.pow(m,ad)) * (1-0.8*cloudiness) * islit)
    else:
        epw = open(epwpath.strip('\n'), 'r')
        for line in epw:
            checklist = line.split(',')
            if checklist[1] == str(month) and checklist[2] == str(day) and len(checklist) == 35:
                # len(checklist) is to make sure that you locate the right line not the file header
                Ebeam_raw.append(float(checklist[14]))
                Ediffuse_raw.append(float(checklist[15]))
        epw.close()
    Ebeam = np.array(Ebeam_raw)
    Ediffuse= np.array(Ediffuse_raw)
    

    
    Ys = []
    for i in range(len(theta_cos)):
        if altitudes[i] < -20:
            Ys.append(0)
        else:
            if theta_cos[i] < -0.2:
                Ys.append(0.45)
            else:
                Ys.append((0.55+0.437*theta_cos[i]+0.313*theta_cos[i]*theta_cos[i]))
    seriesY = np.array(Ys)
    Qdiffuse = Ediffuse * seriesY

    # print('Qdiffuse', Qdiffuse)
    Qbeam = Ebeam * heavisideFilter(theta_cos)
    # print('Qbeam', Qbeam)

    # Calculate load from heat emission
    for i in range(len(seriesTout)):
        Tsols = []
        deltaR = abs(5.67e-8 * math.pow(273.15 + seriesTout[i], 4) - Qbeam[i] - Qdiffuse[i])
        Tsol = seriesTout[i] + 0.15*(Qbeam[i] + Qdiffuse[i]) + deltaR *0.5 /4
        Tsols.append(Tsol)

    seriesTsol = seriesTout
    
    loadSolar = loopSigma(Qbeam + Qdiffuse, RTSsolar)
    loadWall = loopSigma(setTemp - seriesTsol, CTSFwall) * Uwall
    # Calculate load from internal source
    fracRad = 0.6
    Qconv = (ampPeople*schPeople + ampLight*schLight) * (1-fracRad) + ampEquip*schEquip
    Qrad = (ampPeople*schPeople + ampLight*schLight) * fracRad
    loadMass = loopSigma(Qconv, RTSmass) + Qrad
    beta_rad = [math.radians(altitudes[i]) if altitudes[i] > 0 else 0.001 for i in range(24)]
    loadStructure = (WWR*Uglazing*(setTemp-seriesTout) + (1-WWR)*loadWall) * H
    distances = []


    for i in range(24):
        if loadStructure[i] > loadMass[i] * Hsill / math.tan(beta_rad[i]):
            D = Hsill / math.tan(beta_rad[i]) + \
            (loadStructure[i] - loadMass[i] * Hsill / math.tan(beta_rad[i])) / \
            loadSolar[i] / SHGC 
        else:
            D = loadStructure[i] / loadMass[i]
        distances.append(D)

    # print('loadStructure', loadStructure)
    # print('loadMass', loadMass)
    # print('loadSolar', loadSolar)
    distances_mix = []
    for i in range(24):
        if loadMass[i] + loadSolar[i] * SHGC * H * WWR / 10 == 0:
            D = 9999
        else:
            D = loadStructure[i] / (loadMass[i] + loadSolar[i] * SHGC * H * WWR / 10)
        distances_mix.append(D)


    Diso = [round(i, 1) for i in heavisideFilter(distances)]
    Dmix = [round(i, 1) for i in heavisideFilter(distances_mix)]
    print(distances_mix)
    if arraySigma(loadMass) + arraySigma(loadSolar * SHGC * H * WWR) == 0:
        Deven = 9999
    else:
        Deven = arraySigma(loadStructure) / (arraySigma(loadMass) + arraySigma(loadSolar * SHGC * H * WWR) / 10)

    Dewin = arraySigma(loadStructure[8:19]) / (arraySigma(loadMass[8:19]) + arraySigma(loadSolar[8:19] * SHGC * H * WWR) / 10)
    # self.exposure_1 = [Hsill / math.tan(beta_rad[i]) for i in range(24)]
    # self.exposure_2 = [(Hsill+H*WWR) / math.tan(beta_rad[i]) for i in range(24)]

    return Diso, Dmix, Deven, Dewin


def bar_base(data_1, data_2, mean_1, mean_2) -> Bar:
    ceiling = 10
    data_1 = [i if i <=ceiling else ceiling for i in data_1]
    data_2 = [i if i <=ceiling else ceiling for i in data_2]
    c = (
        Bar()
        .add_xaxis([str(i) for i in range(24)])
        .add_yaxis("D_iso", data_1, 
            gap='-100%', color="#fa897b",
            label_opts=opts.LabelOpts(is_show=False),
            itemstyle_opts=opts.ItemStyleOpts(opacity=0.2)
            )
        .add_yaxis("D_mix", data_2, 
            gap='-100%', color="#ccc",
            label_opts=opts.LabelOpts(is_show=False)
            )
        # .add_yaxis("商家B", [randrange(0, 100) for _ in range(24)])
        .set_global_opts(
            title_opts=opts.TitleOpts(),
            tooltip_opts=opts.ToolboxOpts(is_show=False),
            yaxis_opts=opts.AxisOpts(
                name="Offset",
                type_="value",
                min_=0,
                max_=ceiling,
                interval=1,
                # axislabel_opts=opts.LabelOpts(formatter="{value}m"),
                axistick_opts=opts.AxisTickOpts(is_show=True),
                splitline_opts=opts.SplitLineOpts(is_show=True),
            ),
            # xaxis_opts=opts.AxisOpts(name="Hour"),
            legend_opts=opts.LegendOpts(pos_right=50)
        )
        .set_series_opts(
            label_opts=opts.LabelOpts(is_show=False),
            markline_opts=opts.MarkLineOpts(
                is_silent= True,
                symbol= None,
                symbol_size=0.2,
                data=[opts.MarkLineItem(y=mean_1, name="D_mean"), opts.MarkLineItem(y=mean_2, name="D_ewin", coord=[9,mean_2])],
                linestyle_opts=opts.LineStyleOpts(type_='solid', color='#777', opacity=0.6)
            )
        )
    )

    return c


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/testRoute", methods=['GET', 'POST'])
def test():
    if request.method == "POST":
        city = request.json["inp_city"]
        date = request.json["inp_date"]
        psi = float(request.json["inp_axis"])
        try:
            cloudiness = float(request.json["inp_clo"])
        except:
            cloudiness = None
        try:
            maxTemp = float(request.json["inp_temp1"])
        except:
            maxTemp = None
        try:
            minTemp = float(request.json["inp_temp2"])
        except:
            minTemp = None

        # params validated by js
        Uwall = float(request.json["inp_uwall"])
        Uglazing = float(request.json["inp_uwin"])
        ampPeople = float(request.json["inp_peo"])
        ampLight = float(request.json["inp_lgt"])
        ampEquip = float(request.json["inp_eqp"])
        setTemp = float(request.json["inp_spt"])
        SHGC = float(request.json["inp_shgc"])
        ctsf = request.json["inp_ctsf"]
        rts = request.json["inp_rtf"]
        H = float(request.json["inp_h"])
        Hsill = float(request.json["inp_sill"])
        WWR = float(request.json["inp_wwr"])
        
        # optional inputs
        cond = float(request.json["inp_cond"])
        den = float(request.json["inp_den"])
        cp = float(request.json["inp_cp"])
        thickness = float(request.json["inp_thickness"])


        Diso, Dmix, Deven, Dewin = offsetCalc(city, date, psi, cloudiness, maxTemp, minTemp, Uwall, Uglazing, \
            ampPeople, ampLight, ampEquip, setTemp, SHGC, ctsf, rts, H, Hsill, WWR, cond, den, cp, thickness)
        c = bar_base(Diso, Dmix, Deven, Dewin)
        return c.dump_options_with_quotes()