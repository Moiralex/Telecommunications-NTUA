import numpy as np
import matplotlib.pyplot as plt

import math
from sympy.combinatorics.graycode import GrayCode
import statistics as st

from scipy import special

from scipy.io import wavfile
import scipy.io

##################################################################################################################################################
#ΕΡΩΤΗΜΑ 1

#μερος α
def sampling1_2_4(AM,fm,title,fs1): #συναρτηση δειγματοληψίας
	plt.clf()
	
	T=1/(2*fm) #περίοδος σήματος y
	t=np.arange(0,4*T,(1/fs1))
	y=np.cos(2*np.pi*fm*t)*np.cos(2*np.pi*(AM+2)*fm*t)

	plt.figure(figsize=(20,8))
	plt.xlabel("time(s)")
	plt.ylabel("voltage(V)")
	plt.title(title)
	plt.stem(t,y,'b')
	plt.grid()
	plt.savefig(title)
	return
sampling1_2_4(9,4,"fs1=20*fm_el18049",20*4) #εκτέλεση δειγματοληψιας για AM=9,fm=4,fs=20*fm
sampling1_2_4(1,9,"fs1=20*fm_el18081",20*9) #εκτέλεση AM=1,fm=9,fs=20*fm
sampling1_2_4(9,4,"fs2=100*fm_el18049",100*4) #εκτέλεση AM=9,fm=4,fs=100*fm
sampling1_2_4(1,9,"fs2=100*fm_el18081",100*9) #εκτέλεση AM=1,fm=9,fs=100*fm

def sampling3(AM,fm,title,fs1,fs2): #συναρτηση δειγματοληψίας(κοινό διάγραμμα)
	plt.clf()
	T=1/(2*fm)
	t2=np.arange(0,4*T,(1/fs2))
	y2=np.cos(2*np.pi*fm*t2)*np.cos(2*np.pi*(AM+2)*fm*t2)
	plt.figure(figsize=(20,8))
	plt.stem(t2,y2,'b', use_line_collection=True)
	t1=np.arange(0,4*T,(1/fs1))
	y1=np.cos(2*np.pi*fm*t1)*np.cos(2*np.pi*(AM+2)*fm*t1)
	plt.stem(t1,y1,'r', use_line_collection=True)
	plt.xlabel("time(s)")
	plt.ylabel("voltage(V)")
	plt.title(title)
	plt.grid()
	plt.savefig(title)
sampling3(9,4,"both_el18049",20*4, 100*4)
sampling3(1,9,"both_el18081",20*9, 100*9)

#μερος β

sampling1_2_4(9,4,"fs=5*fm_el18049",5*4) #εκτέλεση δειγματοληψιας για AM=9,fm=4,fs=5*fm
sampling1_2_4(1,9,"fs=5*fm_el18081",5*9) #εκτέλεση δειγματοληψιας για AM=1,fm=9,fs=5*fm


##################################################################################################################################################
#ΕΡΩΤΗΜΑ 2

#μερος α
def midriser_with_gray(AM,fm,title,fs1,bits,title2): #κβαντιση του σηματος και αντιστοιχιση στα επιπεδα Gray
  T=1/(2*fm) 
  t=np.arange(0, 4*T,(1/fs1))
  y=np.cos(2*np.pi*fm*t)*np.cos(2*np.pi*(AM+2)*fm*t)

  L=2**bits
  D=(max(y)-min(y))/L #μεγεθος βηματος κβαντισης
  print("θεωρητικο σq για :", title2 , " =",np.sqrt(D**2/12))
  def midriser(x): #συναρτηση κβαντισης
    if(x>=max(y)):
      return D*(math.floor(x/D)-1/2)
    return D*(math.floor(x/D)+1/2)
  quantization_levels=[]
  i = midriser(min(y)) 
  for j in range(L): #πινακας με επιπεδα κβαντισης
    quantization_levels.append(round(i,6))
    i+=D

  yq=list(map(midriser,y)) #κβαντιση σηματος
  for j in range(len(yq)):
    yq[j]=round(yq[j],6)
  gray_values=GrayCode(bits) 
  gray_values=list(gray_values.generate_gray()) #πινακας με ολες τις τιμες του κωδικα Gray
  
  plt.clf()
  plt.figure()
  plt.xlabel("time(s)")
  plt.ylabel("midriser output")
  plt.title(title)
  plt.yticks(quantization_levels,gray_values)
  plt.step(t,yq)
  plt.grid()
  plt.savefig(title)
  return(yq, y, quantization_levels)
(yq_1, y1_1, ql_1)=midriser_with_gray(9,4,"midriser_gray_el18049",20*4,4,"03118049") #03118049
(yq_2, y1_2, ql_2)=midriser_with_gray(1,9,"midriser_gray_el18081",20*9,5,"03118081") #03118081

#μερος β
def standard_deviation_SNR(yq,y1,bits,title):
  print(title)

  err=y1-yq
  
  s1=st.stdev(err[:10]) #τυπικη αποκλιση σφαλματος κβαντισης για 10 πρωτα δειγματα  
  print("σ10=",s1)

  s2=st.stdev(err[:20]) #τυπικη αποκλιση σφαλματος κβαντισης για 20 πρωτα δειγματα 
  print("σ20=",s2)

  sn=st.stdev(err) #τυπικη αποκλιση σφαλματος κβαντισης για ολα τα δειγματα 
  
  def db(x):
    return 10*math.log10(x)

  sum=0 #υπολογισμος ισχυς σηματος
  for i in range(len(y1)):
    sum = sum + (y1[i])**2 
  #P = sum/len(y1)
  #print("P=",P)
  #SNR3 = P/(sn)**2
  #print("SNR_P=", SNR3)
  #print("SNR3_P(db)=", db(SNR3))
  
  sy=st.stdev(y1)

  SNR1=sy/s1
  SNR1=SNR1**2 #SNR κβαντισης στην πρωτη περιπτωση
  print("SNR10=", SNR1)
  print("SNR10(dB)=", db(SNR1))
  SNR2=sy/s2
  SNR2=SNR2**2 #SNR κβαντισης στην δευτερη περιπτωση
  print("SNR20=", SNR2)
  print("SNR20(dB)=", db(SNR2))

  SNR=sy/sn
  SNR=SNR**2
  #print("SNR=",SNR)
  SNR=db(SNR)
  #print("SNR(dB)=",SNR)

standard_deviation_SNR(yq_1,y1_1,4,"03118049:") #03118049
standard_deviation_SNR(yq_2,y1_2,5,"03118081:") #03118081

#μερος γ
def POLARNRZ(fm,yq1,title,bits,ql): #κωδικοποιηση POLAR NRZ του κβαντισμενου σηματος
  T=1/(2*fm) 
  def polar_NRZ(x,Q): #αντιστοιχιση bit στο καταλληλο επιπεδο τασης της κωδικοποιησης
    b=[]
    for i in range(Q):
      if x[i]=="0":
        b.append((-1)*fm)
      else:
        b.append(fm)
    return b

  def polar_NRZ_bitstream(gray_signal, Q): #μετατροπη σηματος με κωδικοποιηση gray σε κωδικοποιηση POLAR NRZ
      stream=[]
      for i in range(len(gray_signal)):
          stream+=polar_NRZ(gray_signal[i],Q)
  
      return stream

  def q2g(quantized_signal,gray_values,quantize_levels):#μετατροπη κβαντισμενου σηματος σε σημα με κωδικοποιηση gray
      a=[]
      for i in range(len(quantized_signal)): 
          x=quantize_levels.index(quantized_signal[i])
          a.append(gray_values[x])
      return a
  gv=GrayCode(bits)
  gv=list(gv.generate_gray())

  a=q2g(yq1, gv, ql)
  #print(a)
  b=polar_NRZ_bitstream(a, bits)
  #print(b)
  k = int(T/0.001) #πληθος bits ανα περιοδο
  b=b[:k]
  n=np.arange(0, T, T/k)
  #print(n)

  plt.clf()
  plt.figure(1)
  plt.step(n, b)
  plt.xlabel("Time(seconds)")
  plt.ylabel("Amplitude(V)")
  plt.title(title)
  plt.grid()
  plt.savefig(title)

POLARNRZ(4,yq_1,"bitstream_POLARNRZ_el18049",4,ql_1) #03118049
POLARNRZ(9,yq_2,"bitstream_POLARNRZ_el18081",5,ql_2) #03118081

##################################################################################################################################################
#ΕΡΩΤΗΜΑ 3

def add_noise_to_signal(signal, SNR, Eb):	#συναρτησης προσθηκης θορυβου στο σημα
  No = Eb/(10**(SNR/10)) #μετατροπη SNR απο dB
    #print("n0=",No)

  noise=np.random.normal(0,np.sqrt(No/2), size=len(signal))+1j*np.random.normal(0,np.sqrt(No/2), size=len(signal))

  r = noise + signal
    
  return r

def BPAM_constellation(A,title1,title2):
#μερος α
  Tb = 0.5  # διαρκεια bit
  #R = 1/Tb  # ρυθμος μεταδοσης bit
  bits = 46

  bitstream = np.random.randint(2, size=46) #τυχαια ακολουθια 46 bits

  bitstream_original = bitstream.copy()
# διαμορφωση BPAM
  BPAM = []
  BPAM.append(0)
  for i in range(bits):
      bitstream[i] = -A if bitstream[i] == 0 else A
      BPAM.append(bitstream[i])
    
  plt.clf()
  t1=np.arange( 0,23.5, 0.5)
  plt.xlabel("time")
  plt.ylabel("Amplitude(V)")
  plt.title('Signal with B-PAM')
  plt.step(t1,BPAM,"g")
  plt.savefig(title1)

#μερος β

#διαγραμμα αστερισμου του BPAM
  Eb = (A**2) * Tb #μεση τιμη ενεργειας bit
  s0 = -np.sqrt(Eb) #σημειο του διαγραμματος αστερισμου που αντιστοιχει στο ψηφιο 0
  s1 = np.sqrt(Eb)  #σημειο του διαγραμματος αστερισμου που αντιστοιχει στο ψηφιο 1
  x = [s0.real, s1.real] 
  y = [s0.imag, s1.imag] 
  plt.clf()
  #plt.figure(figsize=(10,6))
  plt.grid(True)
  plt.plot(x, y, 'o')
  plt.title('Constellation Diagram for B-PAM')
  plt.xlabel('Real Axis')
  plt.ylabel('Imaginary Axis')
  plt.axhline(y=0, linewidth=0.5, color='black')
  plt.axvline(x=0, linewidth=0.5, color='black')
  bits = ['0', '1']
  for i, bit in enumerate(bits):
      plt.annotate(bit, (x[i] + 0.01, y[i] + 0.003))
  plt.axvline(color="red",label="Decision Threshold")
  plt.legend(loc=0,fontsize="small")
  plt.savefig(title2)
  return BPAM,s1,s0,bitstream_original
(BPAM_1,s1_1,s0_1,bs_1)= BPAM_constellation(4,"03118049_BPAM","03118049_constellation")
(BPAM_2,s1_2,s0_2,bs_2)= BPAM_constellation(9,"03118081_BPAM","03118081_constellation")

#μερος γ
def Signal_with_Noise_Constellation(BPAM,s1,s0,title1,title2,SNR1,A):
  BPAM_50x=[]
  Tb = 0.5  # διαρκεια bit
  Eb = (A**2) * Tb
#χωριζουμε καθε παλμο σε 50 μερη για την καλυτερη προσθηκη θορυβου 
  
  for i in range(46):
    i = i+1
    for j in range(50):
     BPAM_50x.append(BPAM[i])
  t1=np.arange(0,23,0.5/50)

  noisy_5 = add_noise_to_signal(BPAM_50x, SNR1,Eb)#προσθηκη θορυβου στο σημα
  plt.figure(figsize = (20, 8))
  plt.clf()
  plt.grid(True)
  plt.xlabel("time")
  plt.ylabel("Amplitude(V)")
  plt.title(title1)
  plt.plot( t1,noisy_5.real, label='SNR1 = Eb/N0 = 5 dB')
  plt.savefig(title1)

#μερος δ
#διαγραμματα αστερισμου μετα την προσθηκη θορυβου
  plt.clf()
  #plt.figure(figsize = (20, 8))
  plt.grid()
  plt.scatter(np.real(noisy_5)*np.sqrt(Tb),np.imag(noisy_5)*np.sqrt(Tb),label="Received signal")
  plt.scatter(s1,0,label="Trasmitted signal") #σημειο αρχικου διαγραμματος αστερισμου
  plt.scatter(s0,0,label="Trasmitted signal")  #σημειο αρχικου διαγραμματος αστερισμου
  plt.axvline(color="black",label="Decision Threshold")
  plt.legend()
  plt.title(title2)
  plt.xlabel("I")
  plt.ylabel("Q")
  if A == 4:
    plt.ylim(-2.5,2.5)
  else:
    plt.ylim(-5,5)
  plt.title(title2)
  plt.savefig(title2)
Signal_with_Noise_Constellation(BPAM_1,s1_1,s0_1,"03118049_awgn+signal_SNR=5","03118049_ConstellationWithNoise_SNR=5",5,4)
Signal_with_Noise_Constellation(BPAM_2,s1_2,s0_2,"03118081_awgn+signal_SNR=5","03118081_ConstellationWithNoise_SNR=5",5,9)
Signal_with_Noise_Constellation(BPAM_1,s1_1,s0_1,"03118049_awgn+signal_SNR=15","03118049_ConstellationWithNoise_SNR=15",15,4)
Signal_with_Noise_Constellation(BPAM_2,s1_2,s0_2,"03118081_awgn+signal_SNR=15","03118081_ConstellationWithNoise_SNR=15",15,9)


#μερος ε
def BER(A,title,s1,s0):
  Tb = 0.5  # διαρκεια bit
  Eb = (A**2) * Tb
  N1=10**6 #αριθμος bits για πειραματικο υπολογισμο
  b=np.random.randint(2,size=N1)

  #πειραματικος υπολογισμος πιθανοτητας λανθασμενου ψηφιου
  c = []
  for i in range(len(b)):
    if b[i] == 0:
      c.append(s0)
    else:
      c.append(s1)
  counter = 0
  probability=[]      
  for d in range(0,16):
    noisy = add_noise_to_signal(c, d, Eb)
    counter = 0
    for i in range(len(noisy)):
       if c[i]*(np.real(noisy[i]))<0:
         counter+=1     
    probability.append(counter/N1)

  #θεωρητικος υπολογισμος πιθανοτητας λανθασμενου ψηφιου
  th=[]
  for d in range(0,16):
	  th.append(0.5*special.erfc(np.sqrt(10**(d/10))))

  plt.clf()
  dB=np.arange(0,16,1)
  plt.yscale("log")
  plt.xlabel("Eb/N0 (dB)")
  plt.ylabel("BER")
  plt.title(title)
  plt.plot(dB, th, '-or', label="Theoretical BPSK")
  plt.plot(dB, probability, 'ob', label="probability BPSK")
  plt.legend()
  plt.grid
  plt.savefig(title)

BER(4,"03118049_BER",s1_1,s0_1)
BER(9,"03118081_BER",s1_2,s0_2)

##################################################################################################################################################
#ΕΡΩΤΗΜΑ 4

#μερος α
 
def QPSK_GRAY_CONSTELLATION(A,title,bs):#παραγωγη διαγραμματος αστερισμου του σηματος με κωδικοποιηση π/4 Gray
  symbol=[]
  sym2plot=[]
  Tb=0.5
  Eb = (A**2)*Tb
  Es = 2*Eb

  for i in range(0, 45, 2):#αντιστοιχιση καθε συμβολο στο καταλληλο συμβολο κωδικοποιησης Gray
    if (bs[i]==0 and bs[i+1]==0):
      symbol.append("00")
      sym2plot.append((np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es))
    elif (bs[i]==0 and bs[i+1]==1):
      symbol.append("01")
      sym2plot.append((np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es))
    elif (bs[i]==1 and bs[i+1]==0):
      symbol.append("10")
      sym2plot.append(-(np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es))
    else:
      symbol.append("11")  
      sym2plot.append(-(np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es))
  
  #τα 4 σημεια στα οποια θα αντιστοιχιζονται τα κωδικοποιημενα κατα π/4 Gray συμβολα
  z00 = (np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)
  z01 = (np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
  z11=-(np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
  z10=-(np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)

  x=[]
  y=[]
  for i in range(23):
    x.append(sym2plot[i].real)
    y.append(sym2plot[i].imag)
  plt.clf()
  plt.figure(figsize=(20,8))
  plt.grid(True)
  plt.plot(x, y, 'o')
  plt.title(title)
  plt.xlabel('Real Axis')
  plt.ylabel('Imaginary Axis')
  plt.axhline(y=0, linewidth=0.5, color='black')
  plt.axvline(x=0, linewidth=0.5, color='black')
  plt.annotate('00', (z00.real+0.15, z00.imag+0.15))
  plt.annotate('01', (z01.real+0.15, z01.imag+0.15))
  plt.annotate('11', (z11.real+0.15, z11.imag+0.15))
  plt.annotate('10', (z10.real+0.15, z10.imag+0.15))
  plt.axhline(color="red")
  plt.axvline(color="red")
  plt.savefig(title)
  return symbol
symbol_1=QPSK_GRAY_CONSTELLATION(4,"03118049_Constellation Diagram for pi4 QPSK",bs_1)
symbol_2=QPSK_GRAY_CONSTELLATION(9,"03118081_Constellation Diagram for pi4 QPSK",bs_2)

#μερος β
#παραγωγη διαγραμματος αστερισμου του σηματος μετα την προσθηκη θορυβο με κωδικοποιηση π/4 Gray
def Constellation_With_Noise(A,SNR1,title,symbol):
  Tb=0.5

  Tb=0.5
  Eb = (A**2)*Tb
  Es = 2*Eb
  z00 = (np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)
  z01 = (np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
  z11=-(np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
  z10=-(np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)

  #αντιστοιχιση bits στα σημεια του διαγραμματος αστερισμου,για καθε bit χρησιμοποιουμε 20 σημεια για καλυτερη απεικονιση του θορυβου
  qpsk_conste_1=[]
  for i in range(23):
    if (symbol[i]=="00"):
      for j in range(20):
        qpsk_conste_1.append(z00)
    elif (symbol[i]=="01"):
      for j in range(20):
        qpsk_conste_1.append(z01)
    elif (symbol[i]=="10"):
      for j in range(20):
        qpsk_conste_1.append(z10)
    else:
      for j in range(20):
        qpsk_conste_1.append(z11)


  plt.clf()
  plt.figure(figsize = (20, 8))
  qpsk_conste=add_noise_to_signal(qpsk_conste_1, SNR1 , Es)#προσθηκη θορυβου στα σημεια του διαγραμματος αστερισμου
  plt.grid()
  plt.scatter(np.real(qpsk_conste),np.imag(qpsk_conste),label="Received signal")
  plt.scatter(z00.real,z00.imag,label="Trasmitted signal")
  plt.scatter(z01.real,z01.imag,label="Trasmitted signal")
  plt.scatter(z10.real,z10.imag,label="Trasmitted signal")
  plt.scatter(z11.real,z11.imag,label="Trasmitted signal")
  plt.axvline(color="black",label="Decision Threshold")
  plt.legend()
  plt.xlabel("I")
  plt.ylabel("Q")
  plt.title(title)
  plt.savefig(title)
Constellation_With_Noise(4,5,"03118049_Constellation_With_Noise_SNR=5",symbol_1)
Constellation_With_Noise(9,5,"03118081_Constellation_With_Noise_SNR=5",symbol_2)
Constellation_With_Noise(4,15,"03118049_Constellation_With_Noise_SNR=15",symbol_1)
Constellation_With_Noise(9,15,"03118081_Constellation_With_Noise_SNR=15",symbol_2)

#μερος γ 
def BER_QPSK_BPSK(A,title1,title2):#θεωρητικος και πειραματικος υπολογισμος πιθανοτητας εσφαλμενου ψηφιου για κωδικοποιηση QPSK και BPSK
  Tb=0.5
  Eb = (A**2)*0.5
  Es = 2*Eb
  z00 = (np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)
  z01 = (np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
  z11=-(np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
  z10=-(np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)

  N1=10**5 #πληθος Bits για πειραματικους υπολογισμους
  N2=int(np.floor(N1/2)) #πληθος συμβολων
  bs=np.random.randint(2,size=N1) #παραγωγη τυχαιας ακολουθιας bits

  #ομαδοποιηση Bit ανα δυο
  symbol=[]
  for i in range(0, N1, 2):
    if (bs[i]==0 and bs[i+1]==0):
      symbol.append("00")
    elif (bs[i]==0 and bs[i+1]==1):
      symbol.append("01")
    elif (bs[i]==1 and bs[i+1]==0):
      symbol.append("10")
    else:
      symbol.append("11")  
  #αντιστοιχιση συμβολων στα σημεια του διαγραμματος αστερισμου
  QPSKconste=[]
  for i in range(N2):
    if (symbol[i]=="00"):
      QPSKconste.append(z00)
    elif (symbol[i]=="01"):
      QPSKconste.append(z01)
    elif (symbol[i]=="10"):
      QPSKconste.append(z10)
    else:
      QPSKconste.append(z11)
  probability=[]
  #πειραματικος υπολογισμος πιθανοτητας εσφαλμενου ψηφιου για τιμες SNR 0-15dB με QPSK διαμορφωση
  for d in range(0,16):
    noisy = add_noise_to_signal(QPSKconste, d, Eb)
    counter = 0
    for i in range(len(noisy)):
        if(((symbol[i]=="00") and ((noisy[i].real<0) or (noisy[i].imag>0))) or ((symbol[i]=="01") and ((noisy[i].real<0) or (noisy[i].imag<0))) or ((symbol[i]=="11") and ((noisy[i].real>0) or (noisy[i].imag<0))) or ((symbol[i]=="10") and ((noisy[i].real>0) or (noisy[i].imag>0)))):
         counter+=1 #προσθηκη ενος λαθους στην περιπτωση λανθασμενου συμβολου
        if (((symbol[i]=="00") and ((noisy[i].real<0) and (noisy[i].imag>0))) or ((symbol[i]=="01") and ((noisy[i].real<0) and (noisy[i].imag<0))) or ((symbol[i]=="11") and ((noisy[i].real>0) and (noisy[i].imag<0))) or ((symbol[i]=="10") and ((noisy[i].real>0) and (noisy[i].imag>0)))):
          counter+=1#προσθηκη αλλου ενος λαθους στην περιπτωση δυο λανθασμενων ψηφιων στο ιδιο συμβολο
           
    probability.append(counter/N1)
  
  #θεωρητικος υπολογισμος πιθανοτητας εσφαλμενου ψηφιου για τιμες SNR 0-15dB με QPSK διαμορφωση
  theoretical=[]
  for d in range(0,16):
	  theoretical.append(0.5*special.erfc(np.sqrt(10**(d/10))))

  b=np.random.randint(2,size=N1)
  BPSK = []
  #αντιστοιχιση συμβολων στα σημεια του διαγραμματος αστερισμου 
  for i in range(N1):
    if(b[i]==1):
      BPSK.append(np.sqrt(Eb))
    else:
      BPSK.append(-np.sqrt(Eb))
  #πειραματικος υπολογισμος πιθανοτητας εσφαλμενου ψηφιου για τιμες SNR 0-15dB με BPSK διαμορφωση
  probabilitybpsk=[]
  for d in range(0,16):
    noisy = add_noise_to_signal(BPSK, d, Eb)
    counter = 0
    for i in range(len(noisy)):
      if (((bool(b[i]==0)) & bool(noisy[i] > 0)) | ((bool(b[i]==1) & bool(noisy[i] < 0)))):
         counter+=1     
    probabilitybpsk.append(counter/N1)
  #θεωρητικος υπολογισμος πιθανοτητας εσφαλμενου ψηφιου για τιμες SNR 0-15dB με BPSK διαμορφωση
  theoretical_bpsk=[]
  for d in range(0,16):
	  theoretical_bpsk.append(0.5*special.erfc(np.sqrt(10**(d/10))))


  plt.clf()
  dB=np.arange(0,16,1)
  plt.xlabel("Eb/N0 (dB)")
  plt.ylabel("BER")
  plt.yscale("log")
  plt.plot(dB, theoretical, '-or', label="Theoretical QPSK")
  plt.plot(dB, probability, 'ob', label="probability QPSK")
  plt.legend()
  plt.grid()
  plt.title(title1)
  plt.savefig(title1)

  plt.clf()
  plt.yscale("log")
  plt.xlabel("Eb/N0 (dB)")
  plt.ylabel("BER")
  plt.plot(dB, theoretical, '-or', label="Theoretical QPSK")
  plt.plot(dB, theoretical_bpsk, '-oy', label="Theoretical BPSK")
  plt.plot(dB, probability, 'ob', label="probability QPSK")
  plt.plot(dB, probabilitybpsk, 'og', label="probability BPSK")
  plt.legend()
  plt.grid()
  plt.title(title2)
  plt.savefig(title2)
BER_QPSK_BPSK(4,"03118049_BER_QPSK","03118049_BERcompare")
BER_QPSK_BPSK(9,"03118081_BER_QPSK","03118081_BERcompare")


#μερος δ
def text_QPSK_transmission(SNR1,title1,title2,title3):
  A=1
  Tb=0.5
  Eb = (A**2)*Tb
  Es = 2*Eb

  #μερος i)

  def ascii_to_bits(text, encoding='ascii'): #μετατροπη ascii σε Bits
      bits=bin(int.from_bytes(text.encode(encoding), 'big'))[2:]
      return bits.zfill(8*((len(bits)+7)//8))

  #ανοιγμα και διαβασμα του αρχειου 
  text=open('shannon_odd.txt', 'r')
  tx=text.read()
  text.close()
  binary=ascii_to_bits(tx)
  binary=binary.replace("0b","")#αφαιρεση των "0b" που προκύπτουν απο την μετατροπη καθε χαρακτηρα ascii σε bits
  
  #μερος ii)

  t=np.arange(0,len(binary),8)
  dec=[]
  dec2=[]
  #μετατροπη δυαδικων αριθμων σε δεκαδικους
  for i in range(0, len(binary), 8):
    k=0
    for j in range(8):
       k = k + int(binary[i+j])*(2**(7-j))
    dec.append(k-128)
    dec2.append(k)
  
  #κβαντιση σηματος 
  L=2**8 #αριθμος επιπεδων 
  D=2/L #μεγεθος βηματος
  def midriser(x):
    return D*(np.floor(x/D)+1/2)
  quantizationlevels=[]
  for i in range(L): #πινακας με επιπεδα κβαντισμου
    quantizationlevels.append(round(i,6))
    #i+=D
  y_gray=list(map(midriser,dec)) #σημα μετα την κβαντιση
  gray_values=GrayCode(8)
  gray_values=list(gray_values.generate_gray()) #πινακας με τιμες της κωδικοποιησης Gray

  plt.clf()
  plt.figure(figsize=(30, 5))
  plt.xlabel("time(s)")
  plt.ylabel("midriser output")
  plt.step(t,y_gray)
  plt.grid()
  
  plt.title(title1)
  plt.savefig(title1)
  
  #μερος iii)
  #σημα με κωδικοποιηση Gray
  signalwithGray=[]
  for i in range(len(dec2)):
    a = dec2[i]
    signalwithGray.append(gray_values[a])


  z00 = (np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)
  z01 = (np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
  z11=-(np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
  z10=-(np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)

  bs=''.join(signalwithGray)
  #ομαδοποιηση Bits ανα 2(συμβολα)
  symbol=[]
  for i in range(0, len(bs), 2):
    if (bs[i]=='0' and bs[i+1]=='0'):
      symbol.append("00")
    elif (bs[i]=='0' and bs[i+1]=='1'):
      symbol.append("01")
    elif (bs[i]=='1' and bs[i+1]=='0'):
      symbol.append("10")
    else:
      symbol.append("11")  
  #αντιστοιχιση συμβολων στα σημεια του αστερισμου
  qpsk_conste_1=[]
  for i in range(len(symbol)):
    if (symbol[i]=="00"):
      qpsk_conste_1.append(z00) 
    elif (symbol[i]=="01"):
      qpsk_conste_1.append(z01)
    elif (symbol[i]=="10"):
      qpsk_conste_1.append(z10)   
    else:
      qpsk_conste_1.append(z11)
  
  #μερος iv)
  #προσθηκη θορυβου στα αντιστοιχα σημεια αστερισμου
  plt.figure(figsize = (20, 8))
  qpsk_conste=add_noise_to_signal(qpsk_conste_1, SNR1 , Es)

  #μερος v)
  #διαγραμματα αστερισμου μετα την προσθηκη του θορυβου 
  plt.grid()
  plt.scatter(np.real(qpsk_conste),np.imag(qpsk_conste),label="Received signal")
  plt.scatter(z00.real,z00.imag,label="Trasmitted signal")
  plt.scatter(z01.real,z01.imag,label="Trasmitted signal")
  plt.scatter(z10.real,z10.imag,label="Trasmitted signal")
  plt.scatter(z11.real,z11.imag,label="Trasmitted signal")
  plt.axvline(color="black",label="Decision Threshold")
  plt.legend()
  plt.xlabel("I")
  plt.ylabel("Q")
  plt.savefig(title2)

  #μερος vi)
  #πειραματικος υπολογισμος πιθανοτητας εσφαλμενου ψηφιου με QPSK διαμορφωση
  probability=[]
  counter = 0
  for i in range(len(symbol)):
    if (((symbol[i]=="00") and ((qpsk_conste[i].real<0) or (qpsk_conste[i].imag>0))) or ((symbol[i]=="01") and ((qpsk_conste[i].real<0) or (qpsk_conste[i].imag<0))) or ((symbol[i]=="11") and ((qpsk_conste[i].real>0) or (qpsk_conste[i].imag<0))) or ((symbol[i]=="10") and ((qpsk_conste[i].real>0) or (qpsk_conste[i].imag>0)))):
      counter+=1 #προσθηκη ενος λαθους στην περιπτωση λανθασμενου συμβολου
    if (((symbol[i]=="00") and ((qpsk_conste[i].real<0) and (qpsk_conste[i].imag>0))) or ((symbol[i]=="01") and ((qpsk_conste[i].real<0) and (qpsk_conste[i].imag<0))) or ((symbol[i]=="11") and ((qpsk_conste[i].real>0) and (qpsk_conste[i].imag<0))) or ((symbol[i]=="10") and ((qpsk_conste[i].real>0) and (qpsk_conste[i].imag>0)))):
      counter+=1  #προσθηκη αλλου ενος λαθους στην περιπτωση δυο λανθασμενων ψηφιων στο ιδιο συμβολο  
  probability.append(counter/len(binary))
  print("probability=", counter/len(binary))
  
  #θεωρητικος υπολογισμος πιθανοτητας εσφαλμενου ψηφιου με QPSK διαμορφωση
  theoretical=0.5*special.erfc(np.sqrt(0.5*10**(SNR1/10)))

  print("Theoreticalprobability=", theoretical)

  #μερος vii)
  #αντιστοιχιση σημειων αστερισμου στα καταλληλα συμβολα
  desymbol=[]
  for i in range(0, len(qpsk_conste), 1):
    if (qpsk_conste[i].real>0 and qpsk_conste[i].imag<0):
      desymbol.append('00')
    elif (qpsk_conste[i].real>0 and qpsk_conste[i].imag>0):
      desymbol.append('01')
    elif (qpsk_conste[i].real<0 and qpsk_conste[i].imag<0):
      desymbol.append('10')
    else:
      desymbol.append('11')  

  #ομαδοποιηση bit ανα 8αδες  
  string=[]
  for i in range(0,len(desymbol),4):
    a=''
    for j in range(4):
      a+=desymbol[i+j]
    string.append(a)
  
  reverse_gray2=[]
  for i in range(len(string)): #αντιστροφη Gray
    reverse_gray=gray_values.index(string[i])#ευρεση θεσης στοιχειου στον πινακα gray και αντιστοιχιση στο καταλληλο επιπεδο κβαντισμου(0-255)
    inter=bin(reverse_gray) #μετατροπη επιπεδου κβαντισμου σε δυαδικο και μετατροπης ολων σε 8-bit ακολουθιες
    temp=10-len(str(inter))
    for j in range(temp):
      reverse_gray2.append('0')
    reverse_gray2.append(bin(reverse_gray))

  reverse_gray3=''.join(reverse_gray2)

  decode=reverse_gray3.replace("0b","")

  def bits_to_ascii(bits): #μετατροπη bits σε συμβολα ascii
    debytes=[bits[8*i:8*i+8] for i in range(len(bits)//8)]
    return "".join([chr(int(b,2)) for b in debytes])

  detext=bits_to_ascii(decode)
  
  file1=open(title3,"w")
  file1.write(detext)
text_QPSK_transmission(5,"binary","signal+awgn_constellation_SNR=5","shannon_odd_SNR=5.txt")
text_QPSK_transmission(15,"binary","signal+awgn_constellation_SNR=15","shannon_odd_SNR=15.txt")

##################################################################################################################################################
#ΕΡΩΤΗΜΑ 5

#μερος α
def Sound_QPSK_transmission(SNR1,title1,title2,title3,title4,title5):
  #διαβασμα αρχειου ηχου
  audio=scipy.io.wavfile.read('soundfile1_lab2.wav')
  
  audio2=[]
  audio2=audio[1] #τιμες του ηχητικου σηματος οπως διαβαζονται απο το αρχειο ηχου
  t1=np.arange(0,len(audio2),1)
  plt.figure(figsize=(20,5))
  plt.plot(t1,audio2)
  plt.grid()
  plt.xlabel("time")
  plt.ylabel("Amplitude")
  plt.title(title1)
  plt.savefig(title1)
  
  #μερος β
  #κβαντιση του σηματος 
  L=2**8 #επιπεδα κβαντισης
  max1=max(audio2)-min(audio2)
  D=max1/L #βημα κβαντισης
  
  def midriser(x): #κβαντιστης
    if(x>=10000):
      return D*(np.floor(x/D)-1/2)
    else:
      return D*(np.floor(x/D)+1/2)

  j=midriser(min(audio2))
  quantized_levels=[]
  for i in range(L): #πινακας με επιπεδα κβαντισμου
    quantized_levels.append(np.round(j,6))
    j+=D

  audio_quantized=list(map(midriser,audio2)) #κβαντισμενο σημα
  for i in range(len(audio_quantized)):
	  audio_quantized[i]=np.round(audio_quantized[i],6)
  def quantized_to_gray(qsig,gray_values,quantized_levels): #αντιστοιχιση επιπεδων κβαντισης σε τιμες κωδικοποιησης Gray
    a=[]
    for i in range(len(qsig)):
      x=quantized_levels.index(qsig[i])
      a.append(gray_values[x])
    return a

  plt.clf()
  plt.figure(figsize=(20, 5))
  plt.title("midrised")
  plt.xlabel("time(s)")
  plt.ylabel("midriser output")
  plt.step(t1,audio_quantized)
  plt.grid()
  plt.title(title2)
  plt.savefig(title2)

  gray_values=GrayCode(8)
  gray_values=list(gray_values.generate_gray()) #πινακας με τιμες της κωδικοποιησης Gray
  
  #κωδικοποιηση Gray του κβαντισμενου σηματος
  audiogray=quantized_to_gray(audio_quantized,gray_values,quantized_levels)
  Es = (1**2)
  Eb = Es/2
  Tb=0.5
  
  #τα 4 σημεια στα οποια θα αντιστοιχιζονται τα κωδικοποιημενα κατα π/4 Gray συμβολα
  z00 = (np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)
  z01 = (np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
  z11=-(np.sqrt(2)/2)*np.sqrt(Es)+1j*(np.sqrt(2)/2)*np.sqrt(Es)
  z10=-(np.sqrt(2)/2)*np.sqrt(Es)-1j*(np.sqrt(2)/2)*np.sqrt(Es)

  bs=''.join(audiogray)
  #μερος γ
  #ομαδοποιηση Bits ανα 2(συμβολα)
  symbol=[]
  for i in range(0, len(bs), 2):
    if (bs[i]=='0' and bs[i+1]=='0'):
      symbol.append("00")
    elif (bs[i]=='0' and bs[i+1]=='1'):
      symbol.append("01")
    elif (bs[i]=='1' and bs[i+1]=='0'):
      symbol.append("10")
    else:
      symbol.append("11")  
  #αντιστοιχιση συμβολων στα σημεια του αστερισμου
  qpsk_conste_1=[]
  for i in range(len(symbol)):
    if (symbol[i]=="00"):
      qpsk_conste_1.append(z00)
    elif (symbol[i]=="01"):
      qpsk_conste_1.append(z01)
    elif (symbol[i]=="10"):
      qpsk_conste_1.append(z10)   
    else:
      qpsk_conste_1.append(z11)

  plt.clf()

  #μερος δ
  qpsk_conste=add_noise_to_signal(qpsk_conste_1, SNR1 , Es)#προσθηκη θορυβου στα αντιστοιχα σημεια αστερισμου
  qpsk_conste1=[]
  for i in range(0, len(qpsk_conste), 2):
	  qpsk_conste1.append(qpsk_conste[i])
  plt.grid()
  
  #διαγραμματα αστερισμου μετα την προσθηκη του θορυβου
  plt.scatter(np.real(qpsk_conste),np.imag(qpsk_conste))
  plt.scatter(z00.real,z00.imag,label="Trasmitted signal")
  plt.scatter(z01.real,z01.imag,label="Trasmitted signal")
  plt.scatter(z10.real,z10.imag,label="Trasmitted signal")
  plt.scatter(z11.real,z11.imag,label="Trasmitted signal")
  plt.axhline(y=0, linewidth=0.5, color='black')
  plt.axvline(x=0, linewidth=0.5, color='black')
  plt.axvline(color='red',label="Decision Threshold")
  plt.legend()
  plt.xlabel("I")
  plt.ylabel("Q")
  plt.title(title3)
  plt.savefig(title3)
  
  #μερος ε
  #πειραματικος υπολογισμος πιθανοτητας εσφαλμενου ψηφιου με QPSK διαμορφωση
  counter = 0
  for i in range(len(symbol)):
    if (((symbol[i]=="00") and ((qpsk_conste[i].real<0) or (qpsk_conste[i].imag>0))) or ((symbol[i]=="01") and ((qpsk_conste[i].real<0) or (qpsk_conste[i].imag<0))) or ((symbol[i]=="11") and ((qpsk_conste[i].real>0) or (qpsk_conste[i].imag<0))) or ((symbol[i]=="10") and ((qpsk_conste[i].real>0) or (qpsk_conste[i].imag>0)))):
      counter+=1 #προσθηκη ενος λαθους στην περιπτωση λανθασμενου συμβολου
    if (((symbol[i]=="00") and ((qpsk_conste[i].real<0) and (qpsk_conste[i].imag>0))) or ((symbol[i]=="01") and ((qpsk_conste[i].real<0) and (qpsk_conste[i].imag<0))) or ((symbol[i]=="11") and ((qpsk_conste[i].real>0) and (qpsk_conste[i].imag<0))) or ((symbol[i]=="10") and ((qpsk_conste[i].real>0) and (qpsk_conste[i].imag>0)))):
      counter+=1  #προσθηκη αλλου ενος λαθους στην περιπτωση δυο λανθασμενων ψηφιων στο ιδιο συμβολο
     
  print("propabibility for", title5, " =", 0.5*counter/len(symbol))

  #θεωρητικος υπολογισμος πιθανοτητας εσφαλμενου ψηφιου με QPSK διαμορφωση
  theoretical=0.5*special.erfc(np.sqrt(0.5*10**(SNR1/10)))

  print("theoretical_probability for", title5, " =", theoretical)
  
  #μερος ζ
  #αντιστοιχιση σημειων αστερισμου στα καταλληλα συμβολα
  desymbol=[]
  for i in range(0, len(qpsk_conste), 1):
    if (qpsk_conste[i].real>0 and qpsk_conste[i].imag<0):
      desymbol.append('00')
    elif (qpsk_conste[i].real>0 and qpsk_conste[i].imag>0):
      desymbol.append('01')
    elif (qpsk_conste[i].real<0 and qpsk_conste[i].imag<0):
      desymbol.append('10')
    else:
      desymbol.append('11') 
  
  #ομαδοποιηση bit ανα 8αδες  
  string=[]
  for i in range(0,len(desymbol),4):
    a=''
    for j in range(4):
      a+=desymbol[i+j]
    string.append(a)

  reverse_gray2=[]
  for i in range(len(string)):
    reverse_gray=gray_values.index(string[i])#ευρεση θεσης στοιχειου στον πινακα gray και αντιστοιχιση στο καταλληλο επιπεδο κβαντισης
    reverse_gray2.append(quantized_levels[reverse_gray])

  reverse_gray2=np.array(reverse_gray2, dtype=np.int16)
  scipy.io.wavfile.write(title4, 44100, reverse_gray2.astype(np.int16))#εξαγωγη αρχειου ηχου

Sound_QPSK_transmission(4,"audio","quantized audio","awgn+audio_constellation_SNR=4",'soundfile_result_SNR=4.wav',"SNR:4")
Sound_QPSK_transmission(14,"audio","quantized audio","awgn+audio_constellation_SNR=14",'soundfile_result_SNR=14.wav',"SNR:14")