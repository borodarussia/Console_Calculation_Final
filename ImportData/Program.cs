using System;
using System.Collections.Generic;
using System.IO;

namespace ImportData
{
    class Program
    {
        public const int NumParameters = 55;
        public const double k = 1.4;
        public const double R = 287;
        public const string FilePath = @"D:\Reposit\C#\Console_Calculation_Final\ImportData\ImportData_iter_1.csv";
        static void Main(string[] args)
        {
            string[] rowImport = File.ReadAllLines(FilePath);
            int rowNumber = rowImport.Length;
            int columnNumber = (rowImport[0].Split(",")).Length;
            string[,] tableString = new string[rowNumber, columnNumber];
            for (int i = 0; i < rowNumber; i++)
            {
                string[] boolRow = rowImport[i].Split(",");
                for (int j = 0; j < boolRow.Length; j++)
                {                    
                    if (boolRow[j] == "")
                    {
                        tableString[i, j] = "0";
                    }
                    else
                    {
                        tableString[i, j] = boolRow[j];
                    }
                }
            }
            double[,] tableDouble = new double[rowNumber, columnNumber - 1];

            for (int i = 0; i < rowNumber; i++)
            {
                for (int j = 1; j < columnNumber; j++)
                {
                    tableDouble[i, j - 1] = double.Parse(tableString[i, j]);
                }
            }

            int NumChannel = Convert.ToInt32(tableDouble[1, 0]);
            int NumSection = NumChannel + 1;
            int NumSectionIn = Convert.ToInt32(tableDouble[8, 0]);
            int NumSectionOut = Convert.ToInt32(tableDouble[9, 0]);
            double[] SectionIn = new double[NumSectionIn];
            double[] SectionOut = new double[NumSectionOut];
            if (NumSectionIn + NumSectionOut <= NumSection)
            {
                for (int i = 0; i < NumSectionIn; i++)
                {
                    SectionIn[i] = Convert.ToInt32(tableDouble[2, i]);
                }
                for (int i = 0; i < NumSectionOut; i++)
                {
                    SectionOut[i] = Convert.ToInt32(tableDouble[2, i]);
                }
            }
            else
            {
                Console.WriteLine("\n\n***\n\nЗаданное количество входных и выходных сечений больше общего количества сечений в системе\n\n***\n\n");
                return;
            }

            double[,] ChannelParameters = new double[NumParameters, NumChannel]; //Объявление массива с параметрами для канала типа трубы
            for (int i = 0; i < NumChannel; i++)
            {
                ChannelParameters[0, i] = Convert.ToInt32(tableDouble[2, i]); // Входное сечение
                ChannelParameters[1, i] = Convert.ToInt32(tableDouble[3, i]); //Выходное сечение
                ChannelParameters[2, i] = tableDouble[4, i] / 1000; // Диаметр на входе
                ChannelParameters[3, i] = tableDouble[5, i] / 1000; // Диаметр на выходе
                ChannelParameters[43, i] = tableDouble[15, i]; // Гидравлическое сопротивление
                /*ChannelParameters[8, i] = Convert.ToDouble(Console.ReadLine()) / 1000; // Длина канала
                ChannelParameters[4, i] = (ChannelParameters[2, i] + ChannelParameters[3, i]) / 2; // Диаметр по средине канала
                ChannelParameters[5, i] = Math.PI * Math.Pow(ChannelParameters[2, i], 2) / 4; // Площадь на входе
                ChannelParameters[6, i] = Math.PI * Math.Pow(ChannelParameters[3, i], 2) / 4; // Площадь на выходе
                ChannelParameters[7, i] = Math.PI * Math.Pow(ChannelParameters[4, i], 2) / 4; // Площадь по средине
                ChannelParameters[48, i] = Functions.DiamG(ChannelParameters[7, i]); // Гидравлический диаметр*/
            }

            /*Console.Write("\n\n***\n\nПересчитывать гидравлическое сопротивление в каналах?(+ Да/- Нет): ");
            string ReCalcKsi = Console.ReadLine();*/

            for (int i = 0; i < NumSectionIn; i++)
            {
                for (int j = 0; j < NumChannel; j++)
                {
                    if (SectionIn[i] == ChannelParameters[0, j])
                    {
                        ChannelParameters[9, j] = tableDouble[12, i];
                        ChannelParameters[12, j] = tableDouble[13, i];
                    }
                }
            }

            for (int i = 0; i < NumSectionOut; i++)
            {
                for (int j = 0; j < NumChannel; j++)
                {
                    if (SectionIn[i] == ChannelParameters[0, j])
                    {
                        ChannelParameters[10, j] = tableDouble[14, i];
                    }
                }
            }


            for (int i = 0; i < NumSectionIn; i++)
            {
                Console.WriteLine(SectionIn[i]);
            }
            for (int i = 0; i < NumSectionOut; i++)
            {
                Console.WriteLine(SectionOut[i]);
            }



            Console.ReadLine();
        }
    }
}
