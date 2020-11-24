# -*- coding: utf-8 -*-
"""
Модуль расчета параметров схемы замещения ВЛ в соответствии с
Руководящими указаниями по релейной защите № 11 Расчеты токов короткого
замыкания для релейной защиты и системной автоматики 110-750 кВ

г.Саратов 24.11.2020г.

"""

import numpy as np
# Применяемые константы
SQRT2=np.sqrt(2)
omegaj = 100j*np.pi # w*j = 2*pi*f*j (f=50 Гц)
#e0 = 8.8541878128e-9 # Электрическая постоянная F/km
C_1_2pie0=17.9751035845223e6 # 1/(2*np.pi*e0)
C_2pij_100=0.062831853071796j # 2*pi/100


class opora:
    '''Класс для хранения информации о высоковольтной опоре воздушной ЛЭП '''
    def __init__(self,name,C1,C2,T1,T2):
        '''Конструктор объекта хранения информации об высоковольтной опоре воздушной ЛЭП
        C1 - Координаты подвеса проводов 1-ой цепи в комплексном виде, м
        C2 - Координаты подвеса проводов 2-ой цепи в комплексном виде, м
        T1 - Координаты подвеса троса №1 в комплексном виде, м
        T2 - Координаты подвеса троса №2 в комплексном виде, м
        Для опор у которых отсутствует 2-ая цепь или 2-ой трос заполнять
        (0j,0j,0j) и 0j соответственно'''
        self.name=name
        self.C1=C1
        self.C2=C2
        self.T1=[T1,]
        self.T2=[T2,]
        self.Parent=None

    def AddHeight(self,name,DY):
        '''Объект опоры созданный с помощью конструктора может выступать
        в качестве родительского объекта для объекта опоры координаты подвеса
        проводов и тросов которого смещены вверх на величину DY, м.
        Результатом метода op.AddHeight(DY) является новый объект класса opora'''
        if self.Parent is None:
            res = opora.__new__(opora)
            res.name=name
            if self.C1[0] != 0:
                res.C1 = [x + DY*1j for x in self.C1]
            else:
                res.C1 = self.C1
            if self.C2[0] != 0:
                res.C2 = [x + DY*1j for x in self.C2]
            else:
                res.C2 = self.C2
            if self.T1[0] != 0:
                res.T1 = [x + DY*1j for x in self.T1]
            else:
                res.T1 = self.T1
            if self.T2[0] != 0:
                res.T2 = [x + DY*1j for x in self.T2]
            else:
                res.T2 = self.T2
            res.Parent = self
            return res
        else:
            return None

    def __repr__(self):
        '''Вывод на экран параметров оборы'''
        return 'Опора {}: \nЦепь 1 - {}; \nЦепь 2 - {}; \nТрос 1 - {}; \nТрос 2 - {}'.format(self.name, self.C1,self.C2,self.T1,self.T2)


class provod:
    '''Класс для хранения информации о марках проводов и тросов'''
    def __init__(self,name,Zud,Diam,Ke=0.95):
        '''Конструктор объекта хранения информации о проводе или тросе
        name - Марка провода
        Zud - Удельное активное сопротивление провода/троса, Ом/км
        Diam - Паспортный диаметр провода/троса, мм
        Ke - Коэффицент эффективного радиуса провода
        Для стальных тросов по рис. 2-21 необходимо задать как активную так и
        индуктивную составляющую сопротивления, последняя используется для
        определения Кэ - коэффицента эффективного радиуса провода'''
        self.name=name
        self.Zud=Zud
        self.Diam=Diam
        if np.imag(Zud) == 0:
            self.Ke = Ke
        else:
            self.Ke = 1/(10**(np.imag(Zud)/0.145))
        self.RadPrEkvC=Diam*0.001/2
        self.RadPrEkvZ=self.Ke*self.RadPrEkvC
        self.Npr = 1
        self.Lp = 0
        self.Parent = None

    def __repr__(self):
        '''Вывод на экран параметров провода/троса'''
        res = 'Провод {}: Zуд = {}; Диаметр = {:5.1f}; Ke = {:5.3f}'.format(self.name, self.Zud,self.Diam,self.Ke)
        if self.Npr>=2:
            res+='; Nпр={}; Lp={}'.format(self.Npr, self.Lp)
        return res

    def RaschPr(self,name,Npr,Lp):
        '''Объект провода/троса созданный с помощью конструктора может выступать
        в качестве родительского для объекта представляющего собой
        расщепленный провод ...
        name - Наименование марки нового расщепленного провода
        Npr - количество проводов
        Lp - длина распорки, м
        Результатом метода pr.RaschPr(Npr,Lp) является новый объект класса provod'''
        if self.Parent is None:
            res = provod.__new__(provod)
            res.name=name
            res.Zud=self.Zud/Npr
            if Npr==2 :
                res.RadPrEkvZ=(self.RadPrEkvZ*Lp)**0.5
                res.RadPrEkvC=(self.RadPrEkvC*Lp)**0.5
            elif Npr==3 :
                res.RadPrEkvZ=(self.RadPrEkvZ*Lp**2)**(1/3)
                res.RadPrEkvC=(self.RadPrEkvC*Lp**2)**(1/3)
            elif Npr==4 :
                res.RadPrEkvZ=(SQRT2*self.RadPrEkvZ*Lp**3)**0.25
                res.RadPrEkvC=(SQRT2*self.RadPrEkvC*Lp**3)**0.25
            res.Diam=2000*res.RadPrEkvC
            res.Ke=self.Ke**(1/Npr)
            res.Npr = Npr
            res.Lp = Lp
            res.Parent = self
            return res
        else:
            return None


class izol:
    '''Класс для хранения информации о марках изоляторов'''
    def __init__(self,name,Len):
        '''Конструктор объекта хранения информации об изоляторе
        name - Марка изолятора
        Len - длина изолятора, м
        Niz - количество изоляторов в гирлянде'''
        self.name=name
        self.Len=Len
        self.Niz=1
        self.LenG=2.5j*Len # (1+1.5)=2.5 длины изолятора
        self.Parent = None

    def __repr__(self):
        '''Вывод на экран параметров изолятора'''
        return 'Изолятор {}: L = {}'.format(self.name, self.Len)

    def __call__(self,Niz):
        '''Формирование нового объекта класса izol на основе исходного с целью
        представления гирлянды изоляторов с числом указанном как аргумент метода
        Применение метода - iz(8), где
        iz - объект класса izol
        8 - количество изоляторов в гирлянде'''
        if self.Parent is None:
            res = izol.__new__(izol)
            res.name='{}x{}'.format(Niz, self.name)
            res.Len=self.Len
            res.Niz=Niz
            res.LenG=1j*res.Len*(Niz+1.5)
            res.Parent = self
            return res
        else:
            return None


class Line:
    '''Класс представляющий трехфазную линию или трос в пределах одного сечения'''

    def __init__(self,sech,name,KOP,OpC,Pr,Iz,RZT=0,fpr=0):#,q1=None,q2=None
        '''Конструктор объекта трехфазной линии или троса в пределах одного сечения
        sech - Ссылка на сечение в которое добавляется трехфазная линия/трос
        name - Наименование трехфазной линии или троса
        KOP - Смещение места установки опоры от начала координат в комплексном виде, м
        OpC - Координаты подвеса проводов трехфазной линии или троса в комплексном виде, м
        подставляется координаты цепи из объекта типа opora
        Pr - Ссылка на тип провода/троса
        Iz - Ссылка на тип изолятора, можно задать непосредственно длину гирлянды, м
        RZT - режим заземления грозозащитного троса:
            0 - трос выступает как отдельная ВЛ
            1 - трос заземлен с одной стороны (влияет только на B0)
            2 - трос заземлен с двух сторон (влияет на Z0 и B0)
        fpr - Стрела провеса провода/троса, м'''
        self.sech = sech# Ссылка на сечение
        self.name = name# Название ВЛ/троса
        self.q1 = None# Узел №1 к которому подключена ВЛ
        self.q2 = None# Узел №2 к которому подключена ВЛ
        self.KOP = KOP# Смещение относительно точки начала отсчета по X и Y в комплекном виде
        self.OpC = OpC# Координаты 3-х точек повеса проводов ВЛ или 1-ой точки повеса троса по конфигурации опоры в комплекном виде
        self.pr_tr=len(OpC)# pr_tr=3 - трехфазная ВЛ, pr_tr=1 - грозозащитный трос
        self.Pr = Pr# Ссылка на Тип провода ВЛ/троса
        self.Iz = Iz# Ссылка на Тип изолятора в гирлянде подвеса провода ВЛ/троса

        #Длина гирлянды изоляторов
        if isinstance(Iz, izol):
            self.Len_girl_izol = Iz.LenG
        else:
            self.Len_girl_izol = Iz*1j
        self.fpr = fpr# Высота провеса проводов в пролете
		#Координаты проводов с учетом смещения опоры и длины гирлянды изоляторов
        dx = KOP - self.Len_girl_izol - (2j/3)*fpr
        self.coord = [x + dx for x in OpC]
        #Расчеты собственных геометрических и электрических параметров трехфазной ВЛ / троса
        if self.pr_tr==3:
            Dp = (np.prod([np.abs(self.coord[x]-self.coord[x-1]) for x in range(3)]))**(1/3)
            SL = (8*np.prod([np.imag(x) for x in self.coord]))**(1/3)
            SM = (np.prod([np.abs(self.coord[x]-np.conj(self.coord[x-1])) for x in range(3)]))**(1/3)
            self.RZT = 0 #Режим заземления троса (Для трехфазной ВЛ RZT = 0)
            zl = self.Pr.Zud.real + self.sech.Rz + C_2pij_100*np.log(self.sech.Dz/self.Pr.RadPrEkvZ)
            zm = self.sech.Rz + C_2pij_100*np.log(self.sech.Dz/Dp)
            al = C_1_2pie0 * np.log(SL/self.Pr.RadPrEkvC)
            am = C_1_2pie0 * np.log(SM/Dp)
        else:
            SL = 2*np.imag(self.coord[0])
            self.RZT = RZT# Режим заземления троса: 0-как отдельная ВЛ, 1-с одной стороны, 2-с двух сторон
            zl = 3 * (self.Pr.Zud.real + self.sech.Rz + C_2pij_100 * np.log(self.sech.Dz/self.Pr.RadPrEkvZ))
            zm = 0
            al = C_1_2pie0 * np.log(SL/self.Pr.RadPrEkvC)
            am = 0
        self.z1 = zl - zm
        self.z0 = zl + 2 * zm
        self.a1 = al - am
        self.a0 = al + 2 * am

        if RZT==0: #Трехфазная ВЛ или трос как отдельная ВЛ
            sech.bp.append(self)
            self.Id=len(sech.bp)
        elif RZT==1: #трос заземлен с одной стороны (влияет только на B0)
            sech.btC.append(self)
        elif RZT==2: #трос заземлен с двух сторон (влияет на Z0 и B0)
            sech.bt.append(self)
            sech.btC.append(self)

    def Mpp(self,pk):
        '''Служебная функция для расчета взаимоиндукции нулевой последовательности
        между проводами двух объектов Line'''
        if self.pr_tr==pk.pr_tr:
            if self.pr_tr==3: #Между двумя трехфазными цепями
                dpp = (np.prod(np.abs([x-y for x in self.coord for y in pk.coord])))**(1/9)
            else: #Между двумя тросами
                dpp = np.abs(self.coord[0]-pk.coord[0])
        else:
            if self.pr_tr==3: #Между трехфазной ВЛ и тросом
                dpp = (np.prod([np.abs(x-pk.coord[0]) for x in self.coord]))**(1/3)
            else: #Между тросом и трехфазной ВЛ
                dpp = (np.prod([np.abs(x-self.coord[0]) for x in pk.coord]))**(1/3)
        return 3*(self.sech.Rz + C_2pij_100*np.log(self.sech.Dz/dpp))

    def App(self,pk):
        '''Служебная функция для расчета взаминого потенциального коэфициента
        нулевой последовательности между проводами двух объектов Line'''
        if self.pr_tr==pk.pr_tr:
            if self.pr_tr==3: #Между двумя трехфазными цепями
                spp = (np.prod(np.abs([x-y for x in self.coord for y in np.conj(pk.coord)])))**(1/9)
                dpp = (np.prod(np.abs([x-y for x in self.coord for y in pk.coord])))**(1/9)
                app = 3 * C_1_2pie0 * np.log(spp/dpp)
            else: #Между двумя тросами
                spp = np.abs(self.coord[0]-np.conj(pk.coord[0]))
                dpp = np.abs(self.coord[0]-pk.coord[0])
                app = C_1_2pie0 * np.log(spp/dpp)
        else:
            if self.pr_tr==3: #Между трехфазной ВЛ и тросом
                spp = (np.prod([np.abs(x-pk.coord[0]) for x in np.conj(self.coord)]))**(1/3)
                dpp = (np.prod([np.abs(x-pk.coord[0]) for x in self.coord]))**(1/3)
                app = C_1_2pie0 * np.log(spp/dpp)
            else: #Между тросом и трехфазной ВЛ
                spp = (np.prod([np.abs(x-self.coord[0]) for x in np.conj(pk.coord)]))**(1/3)
                dpp = (np.prod([np.abs(x-self.coord[0]) for x in pk.coord]))**(1/3)
                app = C_1_2pie0 * np.log(spp/dpp)
        return app

    def ConnectTo(self,q1,q2):
        '''Метод для задания узлов МРТКЗ к которым будет подключаться данная ветвь'''
        self.q1 = q1
        self.q2 = q2

    def __repr__(self):
        '''Вывод на экран параметров трехфазной линии или троса'''
        str1 = 'Линия/трос {} из сечения {}, KOP={}, XYоп = {};'.format(self.name,self.sech.name,self.KOP,self.OpC)
        str2 = '\nПровод {}; Lгирл = {}, РЗТ={}, fпр = {}'.format(self.Pr.name,self.Len_girl_izol,self.RZT,self.fpr)
        return str1+str2

    def __getattr__(self, attrname):
        '''Данный метод позволяет получать результаты расчета по данной линии
        при условии ее режима заземления RZT=0 для линии lk следующим способом:
        для получения Z1 - lk.Z1
        для получения Z0 - lk.Z0
        для получения B1 - lk.B1
        для получения B0 - lk.B0'''
        if self.RZT == 0:
            if attrname in ('Z1', 'z1'):
                attrval = self.sech.Len * self.sech.Z1[self.Id-1,0]
            elif attrname in ('Z0', 'z0'):
                attrval = self.sech.Len * self.sech.Z0[self.Id-1,self.Id-1]
            elif attrname in ('B1', 'b1'):
                attrval = self.sech.Len * self.sech.B1[self.Id-1,0]
            elif attrname in ('B0', 'b0'):
                attrval = self.sech.Len * self.sech.B0[self.Id-1,self.Id-1]
            else:
                return
            return attrval
        return


class sech:
    '''Класс представляющий сечение проводов и тросов для которых может быть
    осуществлен расчет продольных удельных сопротивлений прямой и нулевой последовательностей,
    а также поперечной емкостной проводимости прямой и нулевой последовательностей'''
    def __init__(self,name,Len=1,Rz=0.05,Dz=1000):
        self.name=name
        self.Len=Len
        self.Rz=Rz
        self.Dz=Dz
        self.bp=[]
        self.bt=[]
        self.btC=[]

    def __repr__(self):
        '''Вывод на экран параметров сечения - Названия, длины, Rз и Dз'''
        return 'Сечение {}: L = {}; Rз = {}; Dз = {};'.format(self.name,self.Len,self.Rz,self.Dz)

    def calc(self):
        '''Расчет параметров схемы замещения сечения в соответствии с
        Руководящими указаниями по релейной защите № 11 Расчеты токов короткого
        замыкания для релейной защиты и системной автоматики 110-750 кВ'''
        self.np=len(self.bp)
        self.Z1=np.zeros((self.np,1),dtype=np.cdouble)
        self.Z0=np.zeros((self.np,self.np),dtype=np.cdouble)
        A1=np.zeros((self.np,1),dtype=np.double)
        A0=np.zeros((self.np,self.np),dtype=np.double)
        for ij,pk in enumerate(self.bp):
            self.Z1[ij]=pk.z1
            self.Z0[ij,ij]=pk.z0
            A1[ij]=pk.a1
            A0[ij,ij]=pk.a0
            for ij2,pk2 in enumerate(self.bp[0:ij]):
                M0 = pk.Mpp(pk2)
                self.Z0[ij,ij2]=M0
                self.Z0[ij2,ij]=M0
                aM = pk.App(pk2)
                A0[ij,ij2]=aM
                A0[ij2,ij]=aM
        self.nt=len(self.bt)
        if self.nt>0:
            Ztt=np.zeros((self.nt,self.nt),dtype=np.cdouble)
            Zpt=np.zeros((self.np,self.nt),dtype=np.cdouble)
            Ztp=np.zeros((self.nt,self.np),dtype=np.cdouble)
            for ij,pk in enumerate(self.bt):
                Ztt[ij,ij]=pk.z0
                for ij2,pk2 in enumerate(self.bt[0:ij]):
                    M0 = pk.Mpp(pk2)
                    Ztt[ij,ij2]=M0
                    Ztt[ij2,ij]=M0
                for ij2,pk2 in enumerate(self.bp):
                    M0 = pk.Mpp(pk2)
                    Ztp[ij,ij2]=M0
                    Zpt[ij2,ij]=M0
            self.Z0 -= Zpt @ np.linalg.inv(Ztt) @ Ztp
        self.ntC=len(self.btC)
        if self.ntC>0:
            Att=np.zeros((self.ntC,self.ntC),dtype=np.double)
            Apt=np.zeros((self.np,self.ntC),dtype=np.double)
            Atp=np.zeros((self.ntC,self.np),dtype=np.double)
            for ij,pk in enumerate(self.btC):
                Att[ij,ij]= pk.a0
                for ij2,pk2 in enumerate(self.btC[0:ij]):
                    aM = pk.App(pk2)
                    Att[ij,ij2]=aM
                    Att[ij2,ij]=aM
                for ij2,pk2 in enumerate(self.bp):
                    aM = pk.App(pk2)
                    Atp[ij,ij2]=pk2.pr_tr*aM
                    Apt[ij2,ij]=aM
            A0 -= Apt @ np.linalg.inv(Att) @ Atp
        self.B1 = omegaj/A1
        self.B0 = omegaj*np.linalg.inv(A0)

    def res(self):
        '''Вывод на экран результатов расчета'''
        print(self)
        print('Прямая последовательность')
        print('______________________________________________________________________')
        print('|   |                     |        |        |        |        |      |')
        print('| № |      Название       |   R1   |   X1   |   Z1   |   B1   |  fi  |')
        print('|___|_____________________|________|________|________|________|______|')
        for ij,pk in enumerate(self.bp):
            z = pk.Z1
            R = np.real(z)
            X = np.imag(z)
            Z = np.abs(z)
            fi = 180/np.pi*np.angle(z)
            B = 1e6*np.imag(pk.B1)
            print('|   |                     |        |        |        |        |      |')
            print('|{0:^3}|{1:^21s}|{2:^8.4f}|{3:^8.4f}|{4:^8.4f}|{5:^8.4f}|{6:^6.1f}|'.format(ij+1, pk.name,R,X,Z,B,fi))
        print('|___|_____________________|________|________|________|________|______|')
        print()
        print('Нулевая последовательность')
        print('______________________________________________________________________')
        print('|   |                     |        |        |        |        |      |')
        print('| № |      Название       |   R0   |   X0   |   Z0   |   B0   |  fi  |')
        print('|___|_____________________|________|________|________|________|______|')
        for ij,pk in enumerate(self.bp):
            z = pk.Z0
            R = np.real(z)
            X = np.imag(z)
            Z = np.abs(z)
            fi = 180/np.pi*np.angle(z)
            B = 1e6*np.imag(pk.B0)
            print('|   |                     |        |        |        |        |      |')
            print('|{0:^3}|{1:^21s}|{2:^8.4f}|{3:^8.4f}|{4:^8.4f}|{5:^8.4f}|{6:^6.1f}|'.format(ij+1, pk.name,R,X,Z,B,fi))
            for pk2 in self.bp[ij+1:]:
                ij2 = pk2.Id-1
                z = self.Len*self.Z0[ij,ij2]
                R = np.real(z)
                X = np.imag(z)
                Z = np.abs(z)
                fi = 180/np.pi*np.angle(z)
                B = self.Len*1e6*np.imag(self.B0[ij,ij2])
#                print('|{0:^3s}|{1:^21s}|{2:^8.4f}|{3:^8.4f}|{4:^8.4f}|{5:^8.4f}|{6:^6.1f}|'.format('   ', pk2.name,R,X,Z,B,fi))
                print('|/{0:^2}|{1:^21s}|{2:^8.4f}|{3:^8.4f}|{4:^8.4f}|{5:^8.4f}|{6:^6.1f}|'.format(ij2+1, pk2.name,R,X,Z,B,fi))

        print('|___|_____________________|________|________|________|________|______|')
        print()
