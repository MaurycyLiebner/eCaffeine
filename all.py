import bioread
import matplotlib.pyplot as plt
from scipy import stats
import math
import statistics

def findTWave(data):
    x = data[0]
    y = data[1]

    x0 = int(3*len(y)/5)
    x1 = len(y) - 1
    maxY = max(y[x0:x1])
    maxX = x[y.index(maxY)]
    return [maxX, maxY]

def findRWave(data):
    x = data[0]
    y = data[1]

    maxY = max(y)
    maxX = x[y.index(maxY)]
    return [maxX, maxY]

def findSWave(data):
    x = data[0]
    y = data[1]

    minY = min(y)
    minX = x[y.index(minY)]
    return [minX, minY]

def findQWave(data):
    x = data[0]
    y = data[1]

    r = findRWave(data)
    rx = x.index(r[0])
    prevY = r[1]
    for i in reversed(range(0, rx - 1)):
        yi = y[i]
        if prevY - yi < 0.005 and yi > -0.15 and yi < 0.15:
            return [x[i], yi]
        prevY = yi
    return [0, 0]
    

def findQTInterval(data):
    x = data[0]
    y = data[1]

    q = findQWave(data)
    prevY = q[1]
    for i in reversed(range(0, x.index(q[0]) - 10)):
        yi = y[i]
        if abs(prevY - yi) < 0.0005:
            q = [x[i], yi]
            break
        prevY = yi
    
    t = findTWave(data)
    prevY = t[1]
    for i in range(x.index(t[0]) + 25, len(x)):
        yi = y[i]
        if abs(prevY - yi) < 0.0005 or i == len(x) - 1:
            t = [x[i], yi]
            break
        prevY = yi

    return [q, t]

def calculateAVG(data, iStr = "", title = "", displayWholeSignal = False):
    ecg = data.channels[0]

    if displayWholeSignal:
        plt.figure(iStr + " " + title)
        plt.title(iStr + " " + title)
        plt.xlabel("Czas [s]")
        plt.ylabel("Napięcie [mV]")
        plt.plot(ecg.time_index, ecg.data)

    RIds = []
    RX = []
    RY = []

    w = 8
    ww = 2*w
    x = w
    xMax = len(ecg.data) - ww;

    while x < xMax:
        x0 = x - ww
        x1 = x + ww

        minVal = ecg.data[x0]
        maxVal = ecg.data[x0]
        maxId = x0
        
        for y in range(x0, x1):
            val = ecg.data[y]
            if val > maxVal:
                maxVal = val
                maxId = y
            if val < minVal:
                minVal = val

        firstVal = ecg.data[x0]
        lastVal = ecg.data[x1 - 1]

        if maxVal - minVal > 0.2 and maxId - x0 > w and maxId - x0 < 2*ww - w and maxVal - firstVal > 0.1 and maxVal - lastVal > 0.1:
            if maxId not in RIds:
                RIds.append(maxId)
                RX.append(ecg.time_index[maxId])
                RY.append(ecg.data[maxId])
        
        x += w

    if displayWholeSignal:
        plt.plot(RX, RY, 'o')

    RIds.pop(0)
    RIds.pop(-1)

    width = 0
    timeIndexInc = ecg.time_index[1] - ecg.time_index[0]
    RRDists = []
    for i in range(0, len(RIds) - 1):
        RRDists.append(ecg.time_index[RIds[i + 1]] - ecg.time_index[RIds[i]])
        width += RIds[i + 1] - RIds[i]
    width = 2*(int(round(width/len(RIds)))//2) - 1
    print("\n" + iStr + " " + title + ":")
    print("Width: {0}".format(width))

    timeIndex = []

    # print(timeIndexInc)
    avg = []

    for i in range(0, width):
        avg.append(0)
        timeIndex.append(i*timeIndexInc)

    for i in RIds:
        minId = i - width//2
        maxId = i + width//2
        for j in range(minId, maxId + 1):
            avg[j - minId] += ecg.data[j]/len(RIds)
            
    isoAvg = 0
    nIsoAvg = 20
    for i in range(0, nIsoAvg):
        isoAvg += avg[i]/nIsoAvg

    for i in range(len(avg)):
        avg[i] -= isoAvg
            
    return [timeIndex, avg, statistics.stdev(RRDists)]

def calculateAllAVG(allData, dataId):
    allAvg = []
    allTimeIndex = []
    width = 0
    xyArr = []
    longestX = []
    hr = []
    rrStd = []
    qt = []
    ts = []
    rs = []
    ss = []
    for d in allData:
        data = bioread.read_file(d[dataId])
        xy = calculateAVG(data, str(d[0]))
        width += len(xy[0])
        if len(xy[0]) > len(longestX):
            longestX = xy[0]
        xyArr.append(xy)
        
        hr.append(xy[0][len(xy[0]) - 1])
        rrStd.append(xy[2])
        qtData = findQTInterval(xy)
        qt.append(qtData[1][0] - qtData[0][0])
        ts.append(findTWave(xy)[1])
        rs.append(findRWave(xy)[1])
        ss.append(findSWave(xy)[1])
        
    width = 2*(int(round(width/len(allData)))//2) - 1

    for i in range(0, width):
        allTimeIndex.append(longestX[i])
        allAvg.append(0)
    
    for xy in xyArr:
        for i in range(0, width//2 + 1):
            yId = len(xy[1])//2 - i
            if yId < 0:
                allAvg[len(allAvg)//2 - i] += xy[1][0]/len(allData)                
            else:
                allAvg[len(allAvg)//2 - i] += xy[1][yId]/len(allData)
            if i != 0:
                yId = len(xy[1])//2 + i
                if yId >= len(xy[1]):
                    allAvg[len(allAvg)//2 + i] += xy[1][len(xy[1]) - 1]/len(allData)                
                else:
                    allAvg[len(allAvg)//2 + i] += xy[1][yId]/len(allData)
    plt.figure("Porównanie Wszystkich EKG", figsize=(8, 4))
    plt.title("Porównanie Wszystkich EKG")
    plt.xlabel("Czas [s]")
    plt.ylabel("Napięcie [mV]")
    plt.plot(allTimeIndex, allAvg, "k-" if dataId == 1 else "k:")
    plt.legend(['Przed podaniem kofeiny', '30 minut po podaniu kofeiny'])
    if dataId == 2:
        plt.savefig("../wykresy/1_all.svg")
    
    return {"hr" : hr, "rrStd" : rrStd, "qt" : qt, "ts" : ts, "rs" : rs, "ss" : ss};
        

def plotData(data, iStr, title, style):
    avg = calculateAVG(data, iStr, title, False)
    
    plt.figure("Badany " + iStr + " Porównanie EKG", figsize=(8, 4))
    plt.title("Badany " + iStr + " Porównanie EKG")
    plt.xlabel("Czas [s]")
    plt.ylabel("Napięcie [mV]")
    plt.plot(avg[0], avg[1], style)
    plt.legend(['Przed podaniem kofeiny', '30 minut po podaniu kofeiny'])
    
    # t = findTWave(avg)
    # plt.plot(t[0], t[1], 'o')
    
    # q = findQWave(avg)
    # plt.plot(q[0], q[1], 'o')

    # qt = findQTInterval(avg)
    # q = qt[0]
    # t = qt[1]
    # plt.plot(q[0], q[1], 'o')
    # plt.plot(t[0], t[1], 'o')

dataDir = '/home/ailuropoda/Documents/studia/Proj/eCaffeine/'

allData = [[1, dataDir + 'N_M_P_PROJEKT/ID_1_M_EKG_Przed-L05',
            dataDir + 'N_M_P_PROJEKT/ID_1_M_EKG_Po-L05'], # 1
           [3, dataDir + 'N_M_P_PROJEKT/ID_3_K_EKG_Przed-L05',
            dataDir + 'N_M_P_PROJEKT/ID_3_K_EKG_Po-L05'], # 3
           [4, dataDir + 'N_M_P_PROJEKT/ID_4_M_EKG_Przed2-L05',
            dataDir + 'N_M_P_PROJEKT/ID_4_M_EKG_Po-L05'], # 4
           [5, dataDir + 'N_M_P_PROJEKT/ID_5_K_EKG_Przed-L05_Edited',
            dataDir + 'N_M_P_PROJEKT/ID_5_K_EKG_Po-L05'], # 5
           [6, dataDir + 'N_M_P_PROJEKT/ID_6_M_EKG_Przed-L05',
            dataDir + 'N_M_P_PROJEKT/ID_6_M_EKG_Po-L05'], # 6
           [8, dataDir + 'N_M_P_PROJEKT/ID_8_M_EKG_Przed-L05',
            dataDir + 'N_M_P_PROJEKT/ID_8_M_EKG_Po-L05'], # 8
           [9, dataDir + 'N_M_P_PROJEKT/ID_9_M_EKG_Przed-L05',
            dataDir + 'N_M_P_PROJEKT/ID_9_M_EKG_Po-L05'], # 9
           [10, dataDir + 'N_M_P_PROJEKT/ID_10_M_EKG_Przed-L05',
            dataDir + 'N_M_P_PROJEKT/ID_10_M_EKG_Po-L05'], # 10
           [11, dataDir + 'N_M_P_PROJEKT/ID_11_M_EKG_Przed-L05',
            dataDir + 'N_M_P_PROJEKT/ID_11_M_EKG_Po-L05'], # 11
           [12, dataDir + 'N_M_P_PROJEKT/ID_12_K_EKG_Przed-L05',
            dataDir + 'N_M_P_PROJEKT/ID_12_K_EKG_Po-L05'], # 12
           [13, dataDir + 'N_M_P_PROJEKT/ID_13_M_EKG_Przed-L05',
            dataDir + 'N_M_P_PROJEKT/ID_13_M_EKG_Po-L05'], # 13
           [14, dataDir + 'N_M_P_PROJEKT/ID_14_M_EKG_Przed-L05',
            dataDir + 'N_M_P_PROJEKT/ID_14_M_EKG_Po-L05'], # 14
           [15, dataDir + 'N_M_P_PROJEKT/ID_15_K_EKG_Przed-L05',
            dataDir + 'N_M_P_PROJEKT/ID_15_K_EKG_Po-L05'], # 15
           [16, dataDir + 'N_M_P_PROJEKT/ID_16_K_EKG_Przed-L05',
            dataDir + 'N_M_P_PROJEKT/ID_16_K_EKG_Po-L05'] # 16
           ]

for d in allData:
    index = d[0]
    data = bioread.read_file(d[1])
    plotData(data, str(index), "Przed", "k-")
    data = bioread.read_file(d[2])
    plotData(data, str(index), "Po", "k:")
    plt.savefig("../wykresy/1_" + str(index) + ".svg")

data1 = calculateAllAVG(allData, 1)
hr1 = data1["hr"]
rrStd1 = data1["rrStd"]
qt1 = data1["qt"]
ts1 = data1["ts"]
rs1 = data1["rs"]
ss1 = data1["ss"]
data2 = calculateAllAVG(allData, 2)
hr2 = data2["hr"]
rrStd2 = data2["rrStd"]
qt2 = data2["qt"]
ts2 = data2["ts"]
rs2 = data2["rs"]
ss2 = data2["ss"]

def calculateCorrectedQT(hr, qt):
    qtc = []
    for i in range(0, len(hr)):
        qtc.append(1000*qt[i]/math.sqrt(hr[i]))
    return qtc

qtc1 = calculateCorrectedQT(hr1, qt1)
qtc2 = calculateCorrectedQT(hr2, qt2)

def statisticalTest(name1, name2, name3, data1, data2):
    print("\n---------- Testy statystyczne " + name2)
    
    def checkNormalSW(name, data, alpha):
        pv = stats.shapiro(data).pvalue
        print("\n" + name + " Shapiro-Wilk p-value: {0:.3f}".format(pv))
        normal = pv > alphaSW
        if normal:
            print(name + " ma rozkład normalny")
        else:
            print(name + " nie ma rozkładu normalnego")
        return normal

    alphaSW = 0.05
    data1Normal = checkNormalSW("Przed podaniem kofeiny " + name1, data1, alphaSW);
    data2Normal = checkNormalSW("Po podaniu kofeiny " + name1, data2, alphaSW);
    dataDiff = []
    for i in range(0, len(data1)):
        dataDiff.append(data2[i] - data1[i])

    dataDiffNormal = checkNormalSW("Różnica " + name2, dataDiff, alphaSW);

    print()

    alpha = 0.05
    pv = 0
    if (data1Normal and data2Normal) or dataDiffNormal:
        pv = stats.ttest_rel(data1, data2).pvalue
        print(name1 + " Student's t-test p-value: {0:.3f}".format(pv))
    else:
        pv = stats.wilcoxon(data1, data2).pvalue
        print(name1 + " Wilcoxon test p-value: {0:.3f}".format(pv))

    if pv < alpha:
        print("Mediany " + name3 + " przed i po podaniu kofeiny różnią się na poziomie istotności {0:.2f}.".format(alpha))
    else:
        print("Mediany " + name3 + " przed i po podaniu kofeiny nie różnią się na poziomie istotności {0:.2f}.".format(alpha))


statisticalTest("Odstęp R-R", "Odstępu R-R", "Odstępów R-R", hr1, hr2)

plt.figure("Odstęp R-R")
plt.title("Odstęp R-R")
plt.ylabel("Czas [s]")
plt.boxplot([hr1, hr2], tick_labels=["Przed podaniem kofeiny", "30 minut po podaniu kofeiny"],
            medianprops=dict(color='black'))
plt.savefig("../wykresy/2_1_odstep_r-r.svg")



statisticalTest("Odchylenie Standardowe Odstępu R-R", "Odchylenia Standardowego Odstępu R-R", "Odchyleń Standardowych Odstępów R-R", rrStd1, rrStd2)

plt.figure("Odchylenie Standardowe Odstępu R-R")
plt.title("Odchylenie Standardowe Odstępu R-R")
plt.ylabel("Czas [s]")
plt.boxplot([rrStd1, rrStd2], tick_labels=["Przed podaniem kofeiny", "30 minut po podaniu kofeiny"],
            medianprops=dict(color='black'))
plt.savefig("../wykresy/2_2_odch_std_odstep_r-r.svg")


statisticalTest("Odstęp QT", "Odstępu QT", "Odstępów QT", qt1, qt2)

plt.figure("Odstęp QT")
plt.title("Odstęp QT")
plt.ylabel("Czas [s]")
plt.boxplot([qt1, qt2], tick_labels=["Przed podaniem kofeiny", "30 minut po podaniu kofeiny"],
            medianprops=dict(color='black'))
plt.savefig("../wykresy/2_3_odstep_qt.svg")


statisticalTest("Skorygowany Odstęp QT", "Skorygowanego Odstępu QT", "Skorygowanych Odstępów QT", qtc1, qtc2)

plt.figure("Skorygowany Odstęp QT")
plt.title("Skorygowany Odstęp QT")
plt.ylabel("Czas [ms]")
plt.boxplot([qtc1, qtc2], tick_labels=["Przed podaniem kofeiny", "30 minut po podaniu kofeiny"],
            medianprops=dict(color='black'))
plt.savefig("../wykresy/2_4_skorygowany_odstep_qt.svg")


statisticalTest("Amplituda załamka T", "Amplitudy załamka T", "Amplitud załamka T", ts1, ts2)

plt.figure("Amplituda załamka T")
plt.title("Amplituda załamka T")
plt.ylabel("Napięcie [mV]")
plt.boxplot([ts1, ts2], tick_labels=["Przed podaniem kofeiny", "30 minut po podaniu kofeiny"],
            medianprops=dict(color='black'))
plt.savefig("../wykresy/2_5_amplituda_zalamka_t.svg")


statisticalTest("Amplituda załamka R", "Amplitudy załamka R", "Amplitud załamka R", rs1, rs2)

plt.figure("Amplituda załamka R")
plt.title("Amplituda załamka R")
plt.ylabel("Napięcie [mV]")
plt.boxplot([rs1, rs2], tick_labels=["Przed podaniem kofeiny", "30 minut po podaniu kofeiny"],
            medianprops=dict(color='black'))
plt.savefig("../wykresy/2_6_amplituda_zalamka_r.svg")


statisticalTest("Głębokość załamka S", "Głębokości załamka S", "Głębokości załamka S", ss1, ss2)

plt.figure("Głębokość załamka S")
plt.title("Głębokość załamka S")
plt.ylabel("Napięcie [mV]")
plt.boxplot([ss1, ss2], tick_labels=["Przed podaniem kofeiny", "30 minut po podaniu kofeiny"],
            medianprops=dict(color='black'))
plt.savefig("../wykresy/2_7_amplituda_zalamka_s.svg")


plt.show()
