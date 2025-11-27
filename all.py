import bioread
import matplotlib.pyplot as plt
from scipy import stats

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
    for i in range(0, len(RIds) - 1):
        width += RIds[i + 1] - RIds[i]
    width = 2*(int(round(width/len(RIds)))//2) - 1
    print("\n" + iStr + " " + title + ":")
    print("Width: {0}".format(width))

    timeIndex = []
    timeIndexInc = ecg.time_index[1] - ecg.time_index[0]
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
            
    result = [timeIndex, avg]
    return result

def calculateAllAVG(allData, dataId):
    allAvg = []
    allTimeIndex = []
    width = 0
    xyArr = []
    longestX = []
    hr = []
    for d in allData:
        data = bioread.read_file(d[dataId])
        xy = calculateAVG(data, str(d[0]))
        width += len(xy[0])
        if len(xy[0]) > len(longestX):
            longestX = xy[0]
        xyArr.append(xy)
        
        hr.append(xy[0][len(xy[0]) - 1])
        
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
    plt.figure("Porównanie Wszystkich")
    plt.title("Porównanie Wszystkich")
    plt.xlabel("Czas [s]")
    plt.ylabel("Napięcie [mV]")
    plt.plot(allTimeIndex, allAvg)

    return hr;
        

def plotData(data, iStr, title):
    avg = calculateAVG(data, iStr, title, False)
    avgX = avg[0]
    avgY = avg[1]
    plt.figure(iStr + " Porównanie")
    plt.title(iStr + " Porównanie")
    plt.xlabel("Czas [s]")
    plt.ylabel("Napięcie [mV]")
    plt.plot(avgX, avgY)

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
    plotData(data, str(index), "Przed")
    data = bioread.read_file(d[2])
    plotData(data, str(index), "Po")

hr1 = calculateAllAVG(allData, 1)
hr2 = calculateAllAVG(allData, 2)

def statisticalTest(name1, name2, name3, data1, data2):
    def checkNormalSW(name, data, alpha):
        pv = stats.shapiro(data).pvalue
        print("\n" + name + " Shapiro-Wilk p-value: {0:.3f}".format(pv))
        normal = pv < alphaSW
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
plt.boxplot([hr1, hr2], tick_labels=["Przed podaniem kofeiny", "30 minut po podaniu kofeiny"])

plt.show()
