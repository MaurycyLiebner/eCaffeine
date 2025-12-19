import bioread
import matplotlib.pyplot as plt
from scipy import stats
import math
import statistics

plt.rcParams.update({'figure.max_open_warning': 0})

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

    x0 = int(len(y)/2)
    s = [x[x0], y[x0]]
    prevY = s[1]
    for i in range(x0 + 5, len(x)):
        yi = y[i]
        if prevY - yi < 0.005:
            s = [x[i], yi]
            break
        prevY = yi
    return s

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
    
def findQBeginning(data):
    x = data[0]
    y = data[1]
    q = findQWave(data)
    prevY = q[1]
    for i in reversed(range(0, x.index(q[0]) - 10)):
        yi = y[i]
        if abs(prevY - yi) < 0.0025:
            q = [x[i], yi]
            break
        prevY = yi
    return q
    
def findTEnd(data):
    x = data[0]
    y = data[1]
    t = findTWave(data)
    prevY = t[1]
    for i in range(x.index(t[0]) + 25, len(x)):
        yi = y[i]
        if abs(prevY - yi) < 0.0005 or i == len(x) - 1:
            t = [x[i], yi]
            break
        prevY = yi
    return t
        
def findSEnd(data):
    x = data[0]
    y = data[1]
    s = findSWave(data)
    prevY = s[1]
    for i in range(x.index(s[0]) + 10, len(x)):
        yi = y[i]
        if yi - prevY < 0.005:
            s = [x[i], yi]
            break
        prevY = yi
    return s

def findQTInterval(data):
    q = findQBeginning(data)
    t = findTEnd(data)
    return [q, t]

def findQRSInterval(data):
    q = findQBeginning(data)
    s = findSEnd(data)
    return [q, s]

def calculateAVG(data, iStr = "", title = "", displayWholeSignal = False):
    ecg = data.channels[0]

    if displayWholeSignal:
        plt.figure(iStr + " " + title)
        plt.title(iStr + " " + title)
        plt.xlabel("Time [s]")
        plt.ylabel("Voltage [mV]")
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
    # print("\n" + iStr + " " + title + ":")
    # print("Width: {0}".format(width))

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

def calculateAllAVG(allData, dataId, title, filename):
    allAvg = []
    allTimeIndex = []
    width = 0
    xyArr = []
    longestX = []
    hr = []
    rrStd = []
    qt = []
    qrs = []
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
        qrsData = findQRSInterval(xy)
        qrs.append(qrsData[1][0] - qrsData[0][0])
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
    plt.figure(title + " average ECG", figsize=(8, 4))
    plt.title(title + " average ECG")
    plt.xlabel("Time [s]")
    plt.ylabel("Voltage [mV]")
    plt.plot(allTimeIndex, allAvg, "k-" if dataId == 1 else "k:")
    plt.legend(['Before administering caffeine', '30 minutes after administering caffeine'], loc="upper right", prop={'size':8})
    if dataId == 2:
        plt.savefig("../wykresy/" + filename + ".svg")
    
    return {"hr" : hr, "rrStd" : rrStd, "qt" : qt, "qrs" : qrs, "ts" : ts, "rs" : rs, "ss" : ss};
        

def plotData(data, iStr, title, style):
    avg = calculateAVG(data, iStr, title, False)
    
    plt.figure("Subject " + iStr + " ECG comparison", figsize=(8, 4))
    plt.title("Subject " + iStr + " ECG comparison")
    plt.xlabel("Time [s]")
    plt.ylabel("Voltage [mV]")
    plt.plot(avg[0], avg[1], style)
    plt.legend(['Before administering caffeine', '30 minutes after administering caffeine'], loc="upper right", prop={'size':8})
    
    if False:
        t = findTWave(avg)
        plt.plot(t[0], t[1], 'o')
        q = findQWave(avg)
        plt.plot(q[0], q[1], 'o')

    if False:
        qt = findQTInterval(avg)
        q = qt[0]
        t = qt[1]
        plt.plot(q[0], q[1], 'o')
        plt.plot(t[0], t[1], 'o')

    if False:
        qrs = findQRSInterval(avg)
        q = qrs[0]
        s = qrs[1]
        plt.plot(q[0], q[1], 'o')
        plt.plot(s[0], s[1], 'o')

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
            dataDir + 'N_M_P_PROJEKT/ID_16_K_EKG_Po-L05'], # 16
           [17, dataDir + 'N_M_P_PROJEKT/ID_17_M_EKG_Przed-L05',
            dataDir + 'N_M_P_PROJEKT/ID_17_M_EKG_Po-L05'], # 17
           [18, dataDir + 'N_M_P_PROJEKT/ID_18_M_EKG_PrzedR2-L05',
            dataDir + 'N_M_P_PROJEKT/ID_18_M_EKG_PoR2-L05'], # 18
           [19, dataDir + 'N_M_P_PROJEKT/ID_19_K_EKG_Przed-L05',
            dataDir + 'N_M_P_PROJEKT/ID_19_K_EKG_Po-L05'], # 19
           [20, dataDir + 'N_M_P_PROJEKT/ID_20_K_EKG_Przed-L05',
            dataDir + 'N_M_P_PROJEKT/ID_20_K_EKG_Po-L05'] # 20
           ]

def extractData(allData, indexes):
    result = []
    for subject in allData:
        if subject[0] in indexes:
            result.append(subject)
    return result

for d in allData:
    index = d[0]
    data = bioread.read_file(d[1])
    plotData(data, str(index), "Before", "k-")
    data = bioread.read_file(d[2])
    plotData(data, str(index), "After", "k:")
    plt.savefig("../wykresy/1_" + str(index) + ".svg")


def calculateCorrectedQT(hr, qt):
    qtc = []
    for i in range(0, len(hr)):
        qtc.append(1000*qt[i]/math.sqrt(hr[i]))
    return qtc

def statisticalTest(name, data1, data2):
    print("\n---------- Statistical tests " + name)
    
    def checkNormalSW(name, data, alpha):
        pv = stats.shapiro(data).pvalue
        print("\n" + name + " Shapiro-Wilk p-value: {0:.3f}".format(pv))
        normal = pv > alphaSW
        if normal:
            print(name + " is normally distributed")
        else:
            print(name + " isn't normally distributed")
        return normal

    alphaSW = 0.05
    data1Normal = checkNormalSW("Before administering caffeine " + name, data1, alphaSW);
    data2Normal = checkNormalSW("After administering caffeine " + name, data2, alphaSW);
    dataDiff = []
    for i in range(0, len(data1)):
        dataDiff.append(data2[i] - data1[i])

    dataDiffNormal = checkNormalSW("Difference " + name, dataDiff, alphaSW);

    print()

    alpha = 0.05
    pv = 0
    testValueName = ""
    if (data1Normal and data2Normal) or dataDiffNormal:
        testValueName = "Averages"
        pv = stats.ttest_rel(data1, data2).pvalue
        print(name + " Student's t-test p-value: {0:.3f}".format(pv))
    else:
        testValueName = "Medians"
        pv = stats.wilcoxon(data1, data2).pvalue
        print(name + " Wilcoxon test p-value: {0:.3f}".format(pv))

    if pv < alpha:
        print(testValueName + " of " + name + " before and after administering caffeine have a statistically significant difference {0:.2f}.".format(alpha))
    else:
        print(testValueName + " of " + name + " before and after administering caffeine don't have a statistically significant difference {0:.2f}.".format(alpha))


def testSubjects(allData, title, filename):
    print("\n######### Statistical tests " + title)

    data1 = calculateAllAVG(allData, 1, title, filename)
    hr1 = data1["hr"]
    rrStd1 = data1["rrStd"]
    qt1 = data1["qt"]
    qrs1 = data1["qrs"]
    ts1 = data1["ts"]
    rs1 = data1["rs"]
    ss1 = data1["ss"]
    data2 = calculateAllAVG(allData, 2, title, filename)
    hr2 = data2["hr"]
    rrStd2 = data2["rrStd"]
    qt2 = data2["qt"]
    qrs2 = data2["qrs"]
    ts2 = data2["ts"]
    rs2 = data2["rs"]
    ss2 = data2["ss"]


    qtc1 = calculateCorrectedQT(hr1, qt1)
    qtc2 = calculateCorrectedQT(hr2, qt2)


    statisticalTest("R-R interval " + title, hr1, hr2)

    plt.figure("R-R interval " + title)
    plt.title("R-R interval " + title)
    plt.ylabel("Time [s]")
    plt.boxplot([hr1, hr2], tick_labels=["Before administering caffeine", "30 minutes after administering caffeine"],
                medianprops=dict(color='black'))
    plt.savefig("../wykresy/2_1_" + filename + "_r-r_interval.svg")



    statisticalTest("R-R interval standard deviation " + title, rrStd1, rrStd2)

    plt.figure("R-R interval standard deviation " + title)
    plt.title("R-R interval standard deviation " + title)
    plt.ylabel("Time [s]")
    plt.boxplot([rrStd1, rrStd2], tick_labels=["Before administering caffeine", "30 minutes after administering caffeine"],
                medianprops=dict(color='black'))
    plt.savefig("../wykresy/2_2_" + filename + "_r-r_interval_std_dev.svg")


    statisticalTest("QT interval " + title, qt1, qt2)

    plt.figure("QT interval " + title)
    plt.title("QT interval " + title)
    plt.ylabel("Time [s]")
    plt.boxplot([qt1, qt2], tick_labels=["Before administering caffeine", "30 minutes after administering caffeine"],
                medianprops=dict(color='black'))
    plt.savefig("../wykresy/2_3_" + filename + "_qt_interval.svg")


    statisticalTest("Corrected QT interval " + title, qtc1, qtc2)

    plt.figure("Corrected QT interval " + title)
    plt.title("Corrected QT interval " + title)
    plt.ylabel("Time [ms]")
    plt.boxplot([qtc1, qtc2], tick_labels=["Before administering caffeine", "30 minutes after administering caffeine"],
                medianprops=dict(color='black'))
    plt.savefig("../wykresy/2_4_" + filename + "_corrected_qt_interval.svg")


    statisticalTest("T-wave amplitude " + title, ts1, ts2)

    plt.figure("T-wave amplitude " + title)
    plt.title("T-wave amplitude " + title)
    plt.ylabel("Voltage [mV]")
    plt.boxplot([ts1, ts2], tick_labels=["Before administering caffeine", "30 minutes after administering caffeine"],
                medianprops=dict(color='black'))
    plt.savefig("../wykresy/2_5_" + filename + "_wave_amplitude.svg")


    statisticalTest("R-wave amplitude " + title, rs1, rs2)

    plt.figure("R-wave amplitude " + title)
    plt.title("R-wave amplitude " + title)
    plt.ylabel("Voltage [mV]")
    plt.boxplot([rs1, rs2], tick_labels=["Before administering caffeine", "30 minutes after administering caffeine"],
                medianprops=dict(color='black'))
    plt.savefig("../wykresy/2_6_" + filename + "_r_wave_amplitude.svg")


    statisticalTest("S-wave amplitude " + title, ss1, ss2)

    plt.figure("S-wave amplitude " + title)
    plt.title("S-wave amplitude " + title)
    plt.ylabel("Voltage [mV]")
    plt.boxplot([ss1, ss2], tick_labels=["Before administering caffeine", "30 minutes after administering caffeine"],
                medianprops=dict(color='black'))
    plt.savefig("../wykresy/2_7_" + filename + "_s_wave_amplitude.svg")


    statisticalTest("QRS duration " + title, qrs1, qrs2)

    plt.figure("QRS duration " + title)
    plt.title("QRS duration " + title)
    plt.ylabel("Time [s]")
    plt.boxplot([qrs1, qrs2], tick_labels=["Before administering caffeine", "30 minutes after administering caffeine"],
                medianprops=dict(color='black'))
    plt.savefig("../wykresy/2_8_" + filename + "_qrs_duration.svg")

testSubjects(allData, "All subjects", "all_subjects")
womenData = extractData(allData, [3, 5, 12, 15, 16, 19, 20])
testSubjects(womenData, "Women", "women")
menData = extractData(allData, [1, 4, 6, 8, 9, 10, 11, 13, 14, 17, 18])
testSubjects(menData, "Men", "men")
lowIntakeData = extractData(allData, [1, 3, 4, 8, 10, 12, 13])
testSubjects(lowIntakeData, "Low habitual caffeine intake subjects", "low_intake")
highIntakeData = extractData(allData, [2, 5, 6, 7, 9, 11, 14, 15, 16, 17])
testSubjects(highIntakeData, "High habitual caffeine intake subjects", "high_intake")

plt.show()
