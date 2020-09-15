using System;
using System.Security.Cryptography.X509Certificates;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Globalization;
using System.Threading;
using System.Security.Cryptography;
using System.IO;

namespace ConsoleAppSolution
{
    class Functions
    {
        public static double Lambda(double q_lambda, double k = 1.4)
        {
            double _q_lambda_ref = q_lambda;
            double _k = k;
            double _lambda_max = 1;
            double _lambda_min = 0;
            double _lambda_mid = (_lambda_max + _lambda_min) / 2;
            double _q_lambda_cal = Math.Pow((_k + 1) / 2, 1 / (_k - 1)) * _lambda_mid * Math.Pow(1 - (_k - 1) / (_k + 1) * Math.Pow(_lambda_mid, 2), 1 / _k - 1);
            double _delta = Math.Abs(_q_lambda_ref - _q_lambda_cal);
            while (_delta > 0.0001)
            {
                if (_q_lambda_ref < _q_lambda_cal)
                {
                    _lambda_max = _lambda_mid;
                }
                else if (_q_lambda_ref > _q_lambda_cal)
                {
                    _lambda_min = _lambda_mid;
                }
                else
                {
                    break;
                }
                _lambda_mid = (_lambda_max + _lambda_min) / 2;
                //Console.WriteLine(_lambda_mid);
                //Console.ReadLine();
                _q_lambda_cal = Math.Pow((_k + 1) / 2, 1 / (_k - 1)) * _lambda_mid * Math.Pow(1 - (_k - 1) / (_k + 1) * Math.Pow(_lambda_mid, 2), 1 / _k - 1);
                _delta = Math.Abs(_q_lambda_cal - _q_lambda_ref);
            }
            double _lambda = _lambda_mid;
            return _lambda;
        }
        public static double MassFlow(double totpres1, double totpres2, double area_mid, double tottemp_mid, double ksi, double k = 1.4, double R = 287)
        {
            double _totpres1 = totpres1;
            double _totpres2 = totpres2;
            double _area_mid = area_mid;
            double _ksi = ksi;
            double _tottemp_mid = tottemp_mid;
            double _k = k;
            double _R = R;
            double _m = Math.Sqrt(_k / _R * Math.Pow(2 / (_k + 1), (_k + 1) / (_k - 1)));
            double _massflow, _density_mid, _density_mid_tot;
            if (_totpres1 <= _totpres2)
            {
                _massflow = Math.Sqrt((Math.Pow(_totpres1, 2) - Math.Pow(_totpres1 - 1000, 2)) * Math.Pow(_area_mid, 2) / (_ksi * _R * _tottemp_mid));
                _density_mid_tot = _ksi * Math.Pow(_massflow, 2) / (2 * Math.Pow(_area_mid, 2) * (_totpres1 - (_totpres1 - 1000)));

            }
            else
            {
                _massflow = Math.Sqrt((Math.Pow(_totpres1, 2) - Math.Pow(_totpres2, 2)) * Math.Pow(_area_mid, 2) / (_ksi * _R * _tottemp_mid));
                _density_mid_tot = _ksi * Math.Pow(_massflow, 2) / (2 * Math.Pow(_area_mid, 2) * (_totpres1 - _totpres2));
            }
            //Console.WriteLine(_massflow);
            double _tot_pres_mid = _density_mid_tot * _tottemp_mid * _R;
            double _qlambda_mid = _massflow * Math.Sqrt(_tottemp_mid) / (_m * _tot_pres_mid * _area_mid);
            //Console.WriteLine(_qlambda_mid);
            double _lambda_mid = Lambda(_qlambda_mid);
            _density_mid = epsilon_lambda(_lambda_mid) * _density_mid_tot;
            //Console.WriteLine(_lambda_mid);
            double _machnumber_mid = _lambda_mid * (2 / (_k + 1)) / (1 - (_k - 1) / (_k + 1) * Math.Pow(_lambda_mid, 2));
            if (_totpres1 <= _totpres2)
            {
                _massflow = Math.Sqrt(1000 * 2 * _density_mid * Math.Pow(_area_mid, 2) / (_ksi * (1 + 0.25 * Math.Pow(_machnumber_mid, 2))));
            }
            else
            {
                _massflow = Math.Sqrt((_totpres1 - _totpres2) * 2 * _density_mid * Math.Pow(_area_mid, 2) / (_ksi * (1 + 0.25 * Math.Pow(_machnumber_mid, 2))));
            }
            return _massflow;
        }
        public static double Lambda12(double totpres1, double totpres2, double area_mid, double tottemp_mid, double ksi, double k = 1.4, double R = 287)
        {
            double _totpres1 = totpres1;
            double _totpres2 = totpres2;
            double _area_mid = area_mid;
            double _ksi = ksi;
            double _tottemp_mid = tottemp_mid;
            double _k = k;
            double _R = R;
            double _m = Math.Sqrt(_k / _R * Math.Pow(2 / (_k + 1), (_k + 1) / (_k - 1)));
            double _massflow, _density_mid;
            if (_totpres1 <= _totpres2)
            {
                _massflow = Math.Sqrt((Math.Pow(_totpres1, 2) - Math.Pow(_totpres1 - 1000, 2)) * Math.Pow(_area_mid, 2) / (_ksi * _R * _tottemp_mid));
                _density_mid = _ksi * Math.Pow(_massflow, 2) / (2 * Math.Pow(_area_mid, 2) * (_totpres1 - (_totpres1 - 1000)));

            }
            else
            {
                _massflow = Math.Sqrt((Math.Pow(_totpres1, 2) - Math.Pow(_totpres2, 2)) * Math.Pow(_area_mid, 2) / (_ksi * _R * _tottemp_mid));
                _density_mid = _ksi * Math.Pow(_massflow, 2) / (2 * Math.Pow(_area_mid, 2) * (_totpres1 - _totpres2));
            }
            double _tot_pres_mid = _density_mid * _tottemp_mid * _R;
            double _qlambda_mid = _massflow * Math.Sqrt(_tottemp_mid) / (_m * _tot_pres_mid * _area_mid);
            double _lambda_mid = Lambda(_qlambda_mid);
            return _lambda_mid;
        }
        public static double QLambda12(double totpres1, double totpres2, double area_mid, double tottemp_mid, double ksi, double k = 1.4, double R = 287)
        {
            double _totpres1 = totpres1;
            double _totpres2 = totpres2;
            double _area_mid = area_mid;
            double _ksi = ksi;
            double _tottemp_mid = tottemp_mid;
            double _k = k;
            double _R = R;
            double _m = Math.Sqrt(_k / _R * Math.Pow(2 / (_k + 1), (_k + 1) / (_k - 1)));
            double _massflow, _density_mid, _density_mid_tot;
            if (_totpres1 <= _totpres2)
            {
                _massflow = Math.Sqrt((Math.Pow(_totpres1, 2) - Math.Pow(_totpres1 - 1000, 2)) * Math.Pow(_area_mid, 2) / (_ksi * _R * _tottemp_mid));
                _density_mid_tot = _ksi * Math.Pow(_massflow, 2) / (2 * Math.Pow(_area_mid, 2) * (_totpres1 - (_totpres1 - 1000)));

            }
            else
            {
                _massflow = Math.Sqrt((Math.Pow(_totpres1, 2) - Math.Pow(_totpres2, 2)) * Math.Pow(_area_mid, 2) / (_ksi * _R * _tottemp_mid));
                _density_mid_tot = _ksi * Math.Pow(_massflow, 2) / (2 * Math.Pow(_area_mid, 2) * (_totpres1 - _totpres2));
            }
            //Console.WriteLine(_massflow);
            double _tot_pres_mid = _density_mid_tot * _tottemp_mid * _R;
            double _qlambda_mid = _massflow * Math.Sqrt(_tottemp_mid) / (_m * _tot_pres_mid * _area_mid);
            return _qlambda_mid;
        }
        public static double Density(double totpres1, double totpres2, double area_mid, double tottemp_mid, double ksi, double k = 1.4, double R = 287)
        {
            double _totpres1 = totpres1;
            double _totpres2 = totpres2;
            double _area_mid = area_mid;
            double _ksi = ksi;
            double _tottemp_mid = tottemp_mid;
            double _k = k;
            double _R = R;
            double _m = Math.Sqrt(_k / _R * Math.Pow(2 / (_k + 1), (_k + 1) / (_k - 1)));
            double _massflow, _density_mid;
            if (_totpres1 <= _totpres2)
            {
                _massflow = Math.Sqrt((Math.Pow(_totpres1, 2) - Math.Pow(_totpres1 - 1000, 2)) * Math.Pow(_area_mid, 2) / (_ksi * _R * _tottemp_mid));
                _density_mid = _ksi * Math.Pow(_massflow, 2) / (2 * Math.Pow(_area_mid, 2) * (_totpres1 - (_totpres1 - 1000)));
            }
            else
            {
                _massflow = Math.Sqrt((Math.Pow(_totpres1, 2) - Math.Pow(_totpres2, 2)) * Math.Pow(_area_mid, 2) / (_ksi * _R * _tottemp_mid));
                _density_mid = _ksi * Math.Pow(_massflow, 2) / (2 * Math.Pow(_area_mid, 2) * (_totpres1 - _totpres2));
            }
            return _density_mid;
        }
        public static double Viscosity(double p, double t)
        {
            string FilePath = @"D:\Cloud\OneDrive - ssau.ru\Work\AA_Volkov\Gidro_calculate\air.txt";
            if (File.Exists(FilePath)) ;
            var txt = File.ReadAllText(FilePath);
            string[] table = txt.Split("\n");
            string[] rowTable = table[0].Split(" ");
            string[,] RowColTableStr = new string[table.Length, rowTable.Length];
            int numberRows = table.Length;
            string[] numberColumnsArr = table[2].Split(" ");
            int numberColumns = numberColumnsArr.Length;
            double[,] RowColTableDoub = new double[numberRows, numberColumns];
            for (int i = 0; i < numberRows - 1; i++)
            {
                string[] waste_line = table[i].Split(" ");
                for (int j = 0; j < numberColumns; j++)
                {
                    RowColTableStr[i, j] = waste_line[j];
                }
            }
            for (int i = 0; i < numberRows - 1; i++)
            {
                for (int j = 0; j < numberColumns; j++)
                {
                    if (RowColTableStr[i, j] == "")
                    {
                        RowColTableDoub[i, j] = 0;
                    }
                    else
                    {
                        RowColTableDoub[i, j] = Convert.ToDouble(RowColTableStr[i, j]);
                    }
                }
            }
            double t_max, t_min, p_max, p_min;
            int i_max = numberColumns - 1, i_min = 1, j_max = numberRows - 2, j_min = 1;
            p_max = RowColTableDoub[0, numberColumns - 1];
            p_min = RowColTableDoub[0, 1];
            t_min = RowColTableDoub[1, 0];
            t_max = RowColTableDoub[numberRows - 2, 0];
            if (p > p_min && p < p_max)
            {
                for (int i = 1; i < numberColumns - 1; i++)
                {
                    if (p_max - p > RowColTableDoub[0, i] - p_max && RowColTableDoub[0, i] > p)
                    {
                        p_max = RowColTableDoub[0, i];
                        i_max = i;
                    }
                    if (p - p_min > p - RowColTableDoub[0, i] && RowColTableDoub[0, i] < p)
                    {
                        p_min = RowColTableDoub[0, i];
                        i_min = i;
                    }
                }
            }
            else if (p > p_max)
            {
                p_min = RowColTableDoub[0, numberColumns - 2];
            }
            else if (p < p_min)
            {
                p_max = RowColTableDoub[0, 2];
            }
            if (t > t_min && t < t_max)
            {
                for (int j = 1; j < numberRows - 2; j++)
                {
                    if (t_max - t > RowColTableDoub[j, 0] - t_max && RowColTableDoub[j, 0] > t)
                    {
                        t_max = RowColTableDoub[j, 0];
                        j_max = j;
                    }
                    if (t - t_min > t - RowColTableDoub[j, 0] && RowColTableDoub[j, 0] < t)
                    {
                        t_min = RowColTableDoub[j, 0];
                        j_min = j;
                    }
                }
            }
            else if (t > t_max)
            {
                t_min = RowColTableDoub[numberRows - 3, 0];
            }
            else if (t < t_min)
            {
                t_max = RowColTableDoub[2, 0];
            }
            double v_min = (t - t_max) / (t_min - t_max) * (RowColTableDoub[j_min, i_min] - RowColTableDoub[j_max, i_min]) + RowColTableDoub[j_max, i_min];
            double v_max = (t - t_max) / (t_min - t_max) * (RowColTableDoub[j_min, i_max] - RowColTableDoub[j_max, i_max]) + RowColTableDoub[j_max, i_max];
            double _viscosity = (p - p_max) / (p_min - p_max) * (v_min - v_max) + v_max;
            return _viscosity;
        }
        public static double Reynolds(double velocity, double diamg, double density, double viscosity)
        {
            double _reynolds = velocity * diamg * density / viscosity;
            return _reynolds;
        }
        public static double KsiIn(double reynolds, double area_in, double area1, double length, double diamg, double radius_kr = 0, double flow_angle = 0)
        {
            double re = reynolds;
            double area_otn = area_in / area1;
            double r_kr = radius_kr;
            double eta, tau, ksi_an, _ksi_in;
            double diam_g = diamg;
            double l = length;
            double an = flow_angle * Math.PI / 180;
            if (area_otn < 1)
            {
                if (re <= Math.Pow(10, 4))
                {
                    double a0 = 16 - 0.791 * area_otn + 69.93 * Math.Pow(area_otn, 2) - 92.96 * Math.Pow(area_otn, 3);
                    double a1 = -17.2 - 1.05 * area_otn - 123.3 * Math.Pow(area_otn, 2) + 166.5 * Math.Pow(area_otn, 3);
                    double a2 = 7.7 + 2.41 * area_otn + 69.5 * Math.Pow(area_otn, 2) - 69.25 * Math.Pow(area_otn, 3);
                    double a3 = -1.55 - 1.39 * area_otn - 15.2 * Math.Pow(area_otn, 2) + 22.15 * Math.Pow(area_otn, 3);
                    double a4 = 0.118 + 0.201 * area_otn + 1.16 * Math.Pow(area_otn, 2) - 1.8 * Math.Pow(area_otn, 3);
                    _ksi_in = a0 + a1 * Math.Log10(re) + a2 * Math.Pow(Math.Log10(re), 2) + a3 * Math.Pow(Math.Log10(re), 3) + a4 * Math.Pow(Math.Log10(re), 4);
                }
                else
                {
                    if (r_kr / diam_g <= 0.22)
                    {
                        eta = 0.541 - 7.64 * r_kr / diam_g + 39.2 * Math.Pow(r_kr / diam_g, 2) - 69.4 * Math.Pow(r_kr / diam_g, 3);
                    }
                    else
                    {
                        eta = 0.0125;
                    }
                    if (l / diam_g < 2)
                    {
                        tau = 1.5 - 1.53 * l / diam_g + 0.39 * Math.Pow(l / diam_g, 2);
                    }
                    else
                    {
                        tau = 0;
                    }
                    _ksi_in = eta * (1 - area_otn) + tau * Math.Pow(1 - area_otn, 1.5);
                }
            }
            else
            {
                _ksi_in = 0;
            }
            if (an != 0)
            {
                ksi_an = 1 + 0.6 * Math.Cos(an) + 0.4 * Math.Pow(Math.Cos(an), 2);
            }
            else
            {
                ksi_an = 1;
            }
            _ksi_in = _ksi_in * ksi_an;
            return _ksi_in;
        }
        public static double KsiOut(double reynolds, double area_out, double area2, double length, double diamg)
        {
            double re = reynolds;
            double area_otn = area_out / area2;
            double l = length;
            double m, _ksi_out;
            double diam_g = diamg;
            if (area_otn < 1)
            {
                if (l / diam_g < 10)
                {
                    _ksi_out = Math.Pow(1 - area_otn, 2);
                }
                else
                {
                    if (re <= 2300)
                    {
                        _ksi_out = 2 - 2.66 * area_otn + Math.Pow(area_otn, 2);
                    }
                    else
                    {
                        m = 10.15 + 3.853 * Math.Log10(re) + 0.526 * Math.Pow(Math.Log10(re), 2);
                        _ksi_out = Math.Pow(2 * m + 1, 3) * Math.Pow(m + 1, 3) / (4 * Math.Pow(m, 4) * (m + 3) * (2 * m + 3)) - 2 * area_otn * Math.Pow(2 * m + 1, 2) * (m + 1) / (4 * Math.Pow(m, 2) * (m + 2)) + Math.Pow(area_otn, 2);
                    }
                }
            }
            else
            {
                _ksi_out = 0;
            }
            return _ksi_out;
        }
        public static double KsiTr(double reynolds, double area_mid, double length, double diamg)
        {
            double diam_g = diamg;
            double l = length;
            double re = reynolds;
            double lamb_tr;
            if (re < 2300)
            {
                lamb_tr = 64 / re;
            }
            else if (re >= 2300 && re <= 4000)
            {
                lamb_tr = 0.042671 * Math.Log10(re) - 0.1134;
            }
            else if (re > 4000 && re < Math.Pow(10, 5))
            {
                lamb_tr = 0.3164 * Math.Pow(re, -0.25);
            }
            else
            {
                lamb_tr = Math.Pow(1.8 * Math.Log10(re) - 1.64, -2);
            }
            double _ksi_tr = lamb_tr * l / diam_g;
            return _ksi_tr;
        }
        public static double KsiPov(double alpha, double area_mid, double radius, double diamg)
        {
            double diam_g = diamg;
            double alph = alpha * Math.PI / 180;
            double r = radius;
            double k_alph = 1.95 * (0.01 * alph) - 1.21 * Math.Pow(0.01 * alph, 2);
            double k_f = 1;
            double k_r;
            if (r / diam_g < 0.5)
            {
                k_r = 1.2;
            }
            else if (r / diam_g >= 0.5 && r / diam_g <= 1)
            {
                k_r = 0.21 * Math.Pow(r / diam_g, -2.5);
            }
            else
            {
                k_r = 0.21 * Math.Pow(r / diam_g, -0.5);
            }
            double _ksi_pov = k_alph * k_f * k_r;
            return _ksi_pov;
        }
        public static double KsiSum(double ksi_in = 0, double ksi_out = 0, double ksi_tr = 0, double ksi_pov = 0)
        {
            double _ksi_sum = ksi_in + ksi_out + ksi_tr + ksi_pov;
            return _ksi_sum;
        }
        public static double DiamG(double area)
        {
            double diam_g = Math.Sqrt(4 * area / Math.PI);
            return diam_g;
        }
        public static double pi_lambda(double lambda, double k = 1.4)
        {
            double pi_l = Math.Pow(1 - (k - 1) / (k + 1) * Math.Pow(lambda, 2), k / (k - 1));
            return pi_l;
        }
        public static double tau_lambda(double lambda, double k = 1.4)
        {
            double tau = (1 - (k - 1) / (k + 1) * Math.Pow(lambda, 2));
            return tau;
        }
        public static double epsilon_lambda(double lambda, double k = 1.4)
        {
            double epsilon = Math.Pow((1 - (k - 1) / (k + 1) * Math.Pow(lambda, 2)), 1 / (k - 1));
            return epsilon;
        }
        public static double q_lambda(double lambda, double k = 1.4)
        {
            double q_l = Math.Pow((k + 1) / 2, 1 / (k - 1)) * lambda * Math.Pow(1 - (k - 1) / (k + 1) * Math.Pow(lambda, 2), 1 / (k - 1));
            return q_l;
        }
        public static double machnumber(double lambda, double k = 1.4)
        {
            double m = lambda * (2 / (k + 1)) / (1 - (k - 1) / (k + 1) * Math.Pow(lambda, 2));
            return m;
        }
    }

    class Program
    {
        public const int NumParameters = 55;
        public const double k = 1.4;
        public const double R = 287;
        static void Main(string[] args)
        {
            Console.Write("Введите количество каналов в системе: ");
            int NumChannel = Convert.ToInt32(Console.ReadLine()); // Количество каналов в системе
            int NumSection = NumChannel + 1; // Количество сечений в системе
                                             //Далее следует заполнение номеров узлов/сечений для граничных условий
            Console.Write("\n\n***\n\nВведите количество входных сечений в системе: ");
            int NumSectionIn = Convert.ToInt32(Console.ReadLine());
            Console.Write("\n\n***\n\nВведите количество выходных сечений в системе: ");
            int NumSectionOut = Convert.ToInt32(Console.ReadLine());
            double[] SectionIn = new double[NumSectionIn];
            double[] SectionOut = new double[NumSectionOut];
            if (NumSectionIn + NumSectionOut <= NumSection)
            {
                Console.Write("\n\n***\n\n");
                for (int i = 0; i < NumSectionIn; i++)
                {
                    Console.Write("Введите {0}-й номер входного сечения: ", i + 1);
                    SectionIn[i] = Convert.ToInt32(Console.ReadLine());
                }
                Console.Write("\n\n***\n\n");
                for (int i = 0; i < NumSectionOut; i++)
                {
                    Console.Write("Введите {0}-й номер выходного сечения: ", i + 1);
                    SectionOut[i] = Convert.ToInt32(Console.ReadLine());
                }
            }
            else
            {
                Console.Write("\n\n***\n\n");
                Console.WriteLine("Заданное количество входных и выходных сечений больше общего количества сечений в системе");
                return;
            }

            double[,] ChannelParameters = new double[NumParameters, NumChannel]; //Объявление массива с параметрами для канала типа трубы

            Console.Write("\n\n***\n\nЗаполните номера сечений для каналов\n"); //Происходит заполнение первых двух строк массива - номеров сечений в канале
            for (int i = 0; i < NumChannel; i++)
            {
                Console.Write("\nВведите номер входного сечения для канала №{0}: ", (i + 1));
                ChannelParameters[0, i] = Convert.ToInt32(Console.ReadLine());
                Console.Write("Введите номер выходного сечения для канала №{0}: ", (i + 1));
                ChannelParameters[1, i] = Convert.ToInt32(Console.ReadLine());
            }

            Console.WriteLine("\n\n***\n\nЗаполните диаметры и длину канала");
            for (int i = 0; i < NumChannel; i++)
            {
                Console.Write("\nВведите диаметр[мм] на входе в канал №{0}: ", i + 1);
                ChannelParameters[2, i] = Convert.ToDouble(Console.ReadLine()) / 1000; // Диаметр на входе
                Console.Write("Введите диаметр[мм] на входе в канал №{0}: ", i + 1);
                ChannelParameters[3, i] = Convert.ToDouble(Console.ReadLine()) / 1000; // Диаметр на выходе
                Console.Write("Введите длину[мм] канала №{0}: ", i + 1);
                ChannelParameters[8, i] = Convert.ToDouble(Console.ReadLine()) / 1000; // Длина канала
                ChannelParameters[4, i] = (ChannelParameters[2, i] + ChannelParameters[3, i]) / 2; // Диаметр по средине канала
                ChannelParameters[5, i] = Math.PI * Math.Pow(ChannelParameters[2, i], 2) / 4; // Площадь на входе
                ChannelParameters[6, i] = Math.PI * Math.Pow(ChannelParameters[3, i], 2) / 4; // Площадь на выходе
                ChannelParameters[7, i] = Math.PI * Math.Pow(ChannelParameters[4, i], 2) / 4; // Площадь по средине
                ChannelParameters[48, i] = Functions.DiamG(ChannelParameters[7, i]); // Гидравлический диаметр
                Console.Write($"Введите гидравлическое сопротивление для канала {i+1}: ");
                ChannelParameters[43, i] = Convert.ToDouble(Console.ReadLine());
            }

            Console.Write("\n\n***\n\nПересчитывать гидравлическое сопротивление в каналах?(+ Да/- Нет): ");
            string ReCalcKsi = Console.ReadLine(); 

            Console.WriteLine("\n\n***\n\nВведите термогазодинамические параметры для входных сечений");
            for (int i = 0; i < NumSectionIn; i++)
            {
                for (int j = 0; j < NumChannel; j++)
                {
                    if (SectionIn[i] == ChannelParameters[0, j])
                    {
                        Console.Write("\nВведите полное давление для {0}-го сечения: ", SectionIn[i]);
                        ChannelParameters[9, i] = Convert.ToDouble(Console.ReadLine());
                        Console.Write("Введите полную температуру для {0}-го сечения: ", SectionIn[i]);
                        ChannelParameters[12, i] = Convert.ToDouble(Console.ReadLine());
                    }
                }
            }

            Console.WriteLine("\n\n***\n\nВведите термогазодинамические параметры для выходных сечений");
            for (int i = 0; i < NumSectionOut; i++)
            {
                for (int j = 0; j < NumChannel; j++)
                {
                    if (SectionIn[i] == ChannelParameters[0, j])
                    {
                        Console.Write("\nВведите полное давление для {0}-го сечения: ", SectionOut[i]);
                        ChannelParameters[10, i] = Convert.ToDouble(Console.ReadLine());
                        /* Console.Write("Введите полную температуру для {0}-го сечения: ", SectionOut[i]);
                         * ChannelParameters[12, i] = Convert.ToDouble(Console.ReadLine());
                         */

                    }
                }
            }
            //Объявление давлений, необходимых для заполнения промежуточных давлений в каналах
            double pressure_max1 = 0;
            double temper_max = 0;
            double pressure_min1 = Math.Pow(10, 20);

            for (int i = 0; i < NumChannel; i++)
            {
                if (pressure_max1 < ChannelParameters[10, i])
                {
                    pressure_max1 = ChannelParameters[10, i];
                }
                if (pressure_max1 < ChannelParameters[9, i])
                {
                    pressure_max1 = ChannelParameters[9, i];
                }

                if (pressure_min1 > ChannelParameters[10, i] && ChannelParameters[10, i] != 0)
                {
                    pressure_min1 = ChannelParameters[10, i];
                }
                if (pressure_min1 > ChannelParameters[9, i] && ChannelParameters[9, i] != 0)
                {
                    pressure_min1 = ChannelParameters[9, i];
                }

                if (temper_max < ChannelParameters[12, i])
                {
                    temper_max = ChannelParameters[12, i];
                }
            }

            double pressure_max = 0;
            double pressure_min = Math.Pow(10, 20);

            for (int i = 0; i < NumChannel; i++)
            {
                if (NumSectionIn > 1)
                {
                    if (pressure_max < ChannelParameters[10, i] && pressure_max1 != ChannelParameters[10, i])
                    {
                        pressure_max = ChannelParameters[10, i];
                    }
                    if (pressure_max < ChannelParameters[9, i] && pressure_max1 != ChannelParameters[9, i])
                    {
                        pressure_max = ChannelParameters[9, i];
                    }
                }
                else if (NumSectionIn == 1)
                {
                    pressure_max = pressure_max1;
                }

                if (NumSectionOut > 1)
                {
                    if (pressure_min > ChannelParameters[10, i] && ChannelParameters[10, i] != 0 && pressure_min1 != ChannelParameters[10, i])
                    {
                        pressure_min = ChannelParameters[10, i];
                    }
                    if (pressure_min > ChannelParameters[9, i] && ChannelParameters[9, i] != 0 && pressure_min1 != ChannelParameters[9, i])
                    {
                        pressure_min = ChannelParameters[9, i];
                    }
                }
                else if (NumSectionOut == 1)
                {
                    pressure_min = pressure_min1;
                }
            }

            for (int i = 0; i < NumChannel; i++)
            {
                ChannelParameters[13, i] = temper_max;
                ChannelParameters[12, i] = temper_max;
                ChannelParameters[14, i] = temper_max;
                if (ChannelParameters[10, i] == 0)
                {
                    ChannelParameters[10, i] = pressure_min;
                }
                if (ChannelParameters[9, i] == 0)
                {
                    ChannelParameters[9, i] = pressure_min;
                }
            }

            // Определение количества сечений в системе
            int SumSect;
            int[,] SumSection = new int[2, NumSection];

            for (int i = 0; i < NumSection; i++)
            {
                SumSect = 0;
                for (int j = 0; j < NumChannel; j++)
                {
                    if (i + 1 == ChannelParameters[0, j])
                    {
                        SumSect += 1;
                    }
                    else if (i + 1 == ChannelParameters[1, j])
                    {
                        SumSect += 1;
                    }
                }
                SumSection[0, i] = i + 1;
                SumSection[1, i] = SumSect;
            }

            double delta;
            double p_out, p_in;
            double p_max, p_min, p_mid;
            double MassFlowNumerator;

            for (int iter = 0; iter < 100; iter++)
            {
                for (int num_iter = 0; num_iter < 100; num_iter++)
                {
                    p_max = 0;
                    p_min = Math.Pow(10, 20);

                    for (int i = 0; i < NumSection; i++)
                    {
                        delta = 1000000;
                        if (SumSection[1, i] > 1)
                        {
                            for (int j = 0; j < NumChannel; j++)
                            {
                                if (SumSection[0, i] == ChannelParameters[1, j] || SumSection[0, i] == ChannelParameters[0, j])
                                {
                                    if (p_max < ChannelParameters[10, j])
                                    {
                                        p_max = ChannelParameters[10, j];
                                    }
                                    if (p_max < ChannelParameters[9, j])
                                    {
                                        p_max = ChannelParameters[9, j];
                                    }
                                    if (p_min > ChannelParameters[10, j])
                                    {
                                        p_min = ChannelParameters[10, j];
                                    }
                                    if (p_min > ChannelParameters[9, j])
                                    {
                                        p_min = ChannelParameters[9, j];
                                    }
                                }
                            }
                            p_mid = (p_min + p_max) / 2;
                            while (Math.Abs(delta) > 0.0000001)
                            {

                                MassFlowNumerator = 0;
                                for (int j = 0; j < NumChannel; j++)
                                {
                                    if (SumSection[0, i] == ChannelParameters[1, j])
                                    {
                                        p_out = p_mid;
                                        p_in = ChannelParameters[9, j];
                                        ChannelParameters[10, j] = p_out;
                                        MassFlowNumerator += Functions.MassFlow(p_in, p_out, ChannelParameters[7, j], ChannelParameters[14, j], ChannelParameters[43, j]);
                                    }
                                    else if (SumSection[0, i] == ChannelParameters[0, j])
                                    {
                                        p_in = p_mid;
                                        p_out = ChannelParameters[10, j];
                                        ChannelParameters[9, j] = p_in;
                                        MassFlowNumerator -= Functions.MassFlow(p_in, p_out, ChannelParameters[7, j], ChannelParameters[14, j], ChannelParameters[43, j]);
                                    }
                                }
                                if (MassFlowNumerator == delta)
                                {
                                    break;
                                }
                                else
                                {
                                    delta = MassFlowNumerator;
                                    if (delta > 0)
                                    {
                                        p_min = p_mid;
                                        p_mid = (p_max + p_min) / 2;
                                    }
                                    else if (delta < 0)
                                    {
                                        p_max = p_mid;
                                        p_mid = (p_max + p_min) / 2;
                                    }
                                }
                            }

                        }
                    }
                }
                for (int i = 0; i < NumChannel; i++)
                {
                    ChannelParameters[42, i] = Functions.MassFlow(ChannelParameters[9, i], ChannelParameters[10, i], ChannelParameters[7, i], ChannelParameters[14, i], ChannelParameters[43, i]); // Расхода в канале
                    //ChannelParameters[26, i] = Functions.Density(ChannelParameters[9, i], ChannelParameters[10, i], ChannelParameters[7, i], ChannelParameters[14, i], ChannelParameters[43, i]); // Плотность i,j
                    //ChannelParameters[47, i] = ChannelParameters[42, i] / (ChannelParameters[26, i] * ChannelParameters[7, i]); // Скорость для среднего сечения i,j
                    ChannelParameters[29, i] = Functions.Lambda12(ChannelParameters[9, i], ChannelParameters[10, i], ChannelParameters[7, i], ChannelParameters[14, i], ChannelParameters[43, i]); // Лямбда для среднего сечения i,j
                    ChannelParameters[54, i] = Functions.QLambda12(ChannelParameters[9, i], ChannelParameters[10, i], ChannelParameters[7, i], ChannelParameters[14, i], ChannelParameters[43, i]);
                    ChannelParameters[47, i] = ChannelParameters[29, i] * Math.Sqrt(2 * k / (k + 1) * R * ChannelParameters[14, i]); // Скорость для среднего сечения i,j
                    ChannelParameters[26, i] = ChannelParameters[42, i] / (ChannelParameters[47, i] * ChannelParameters[7, i]); // Плотность i,j
                    ChannelParameters[38, i] = Functions.tau_lambda(ChannelParameters[29, i]); // Тау от лямбды i,j
                    ChannelParameters[41, i] = Functions.pi_lambda(ChannelParameters[29, i]); // Пи от лямбды i,j
                    ChannelParameters[23, i] = ChannelParameters[38, i] * ChannelParameters[14, i]; // T i,j
                    ChannelParameters[20, i] = ChannelParameters[26, i] * ChannelParameters[23, i] * R; // статическое давление i,j
                    ChannelParameters[51, i] = Functions.Viscosity(ChannelParameters[20, i], ChannelParameters[23, i]); // Вязкость i, j
                    ChannelParameters[44, i] = Functions.Reynolds(ChannelParameters[47, i], ChannelParameters[48, i], ChannelParameters[26, i], ChannelParameters[51, i]); // Рейнольдс
                    if (ReCalcKsi == "+")
                    {
                        double ksiin = 0, ksiout = 0, ksitr;
                        for (int j = 0; j < NumChannel; j++)
                        {
                            if (ChannelParameters[0, i] == ChannelParameters[1, j])
                            {
                                ksiin += Functions.KsiIn(ChannelParameters[44, i], ChannelParameters[5, i], ChannelParameters[6, j], ChannelParameters[8, i], ChannelParameters[48, i]);
                            }
                            else
                            {
                                ksiin += 0;
                            }
                            if (ChannelParameters[1, i] == ChannelParameters[0, j])
                            {
                                ksiout += Functions.KsiOut(ChannelParameters[44, i], ChannelParameters[6, i], ChannelParameters[5, j], ChannelParameters[8, i], ChannelParameters[48, i]);
                            }
                            else
                            {
                                ksiout += 0;
                            }
                        }
                        ksitr = Functions.KsiTr(ChannelParameters[44, i], ChannelParameters[7, i], ChannelParameters[8, i], ChannelParameters[48, i]);
                        ChannelParameters[43, i] = ksitr + ksiin + ksiout;
                    }                    
                }
            }
            Console.WriteLine("\n\n***\n\n");
            for (int i = 0; i < NumChannel; i++)
            {
                Console.WriteLine($"p1: {ChannelParameters[9, i]} p2: {ChannelParameters[10, i]} G: {ChannelParameters[42, i]} ksi: {ChannelParameters[43, i]} rho: {ChannelParameters[26, i]}");
                Console.WriteLine($"c12: {ChannelParameters[47, i]} lamb12: {ChannelParameters[29, i]} qlamb12: {ChannelParameters[54, i]}");
            }
        }
    }
}
