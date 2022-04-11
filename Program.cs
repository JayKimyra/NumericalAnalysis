using System;
using System.Collections.Generic;
using System.Globalization;

public class NumericalMethods
{

    public delegate double function(double x);
    public delegate double function2(double x, double y);
    public delegate double function3(double x, double y, double z);

    
    public static void Main()
    {
        CultureInfo.CurrentCulture = CultureInfo.GetCultureInfo("en-US");

        function3 F = delegate (double x, double y, double z)
        {
            return -0.001 * (1 / 1.5 - x / (1.5 * 1.5)) * Math.Pow((1 + z * z), 1.5);
        };
        var ans = Runge_Kutta_2nd_Order((x, y, z) => z, F, 0, 1, 500, 0, 0, true);



       
















        var Z = Runge_Kutta((x, y) => -0.001 * (1 / 1.5 - x / (1.5 * 1.5)) * Math.Pow((1 + y * y), 1.5), 0, 1, 500, 0);

        Console.WriteLine("BEGIN");
        var coeff = splineCoeff(Z[0], Z[1]);
        function2 Yx = delegate (double x, double y)
        {
            int n = Z[0].Length;
            for (int i = 0; i < n - 1; i++)
            {
                if (Z[0][i] <= x && Z[0][i + 1] >= x)
                {
                    return coeff[4 * i] + coeff[4 * i + 1] * (x - Z[0][i]) + coeff[4 * i + 2] * (x - Z[0][i]) * (x - Z[0][i]) + coeff[4 * i + 3] * (x - Z[0][i]) * (x - Z[0][i]) * (x - Z[0][i]);
                }
            }
            return 0;
        };
        var ANS = Runge_Kutta((x, y) => Yx(x, y), 0, 1, 500, 0, true);
    }

    //ДиффУры

    static double[] Pikar(function2 f, double a, double b, double h,double y0)
    {
        int n =(int) ((b - a) / h);
        double[] X = new double[n+1];
        for (int i = 0; i < X.Length; i++)
        {
            X[i] = a + i * h;
        }
        double[] Y = new double[n+1];
        Y[0] = y0;
        for (int i = 1; i < Y.Length; i++)
        {
            function fx = x => f(x, Y[i-1]);
            Y[i] = y0 + Integral(fx, a, X[i], n*n,methodName.Middle);
        }
        return Y;
    }

    static double[] Euler(function2 f, double a, double b, double h, double y0)
    {
        int n = (int)((b - a) / h);
        double[] X = new double[n + 1];
        double[] Y = new double[n + 1];
        Y[0] = y0;
        for (int i = 0; i < Y.Length - 1; i++)
        {
            X[i] = a + i * h;
            Y[i+1] = Y[i] + h*f(X[i], Y[i]);
        }
        return Y;
    }

    static double[] Usov_Euler(function2 f, double a, double b, double h, double y0)
    {
        int n = (int)((b - a) / h);
        double[] X = new double[n + 1];
        double[] Y = new double[n + 1];
        Y[0] = y0;
        for (int i = 0; i < Y.Length - 1; i++)
        {
            X[i] = a + i * h;
            Y[i + 1] = Y[i] + h * f(X[i]+h/2, Y[i]+(h/2)* f(X[i], Y[i]));
        }
        return Y;
    }

    static double[][] Runge_Kutta(function2 f, double a, double b, int N, double y0, bool DisplayMode = false)
    {
        double[][] RK = new double[2][];
        double h = (b - a) / N;
        double[] X = new double[N + 1];
        double[] Y = new double[N + 1];
        Y[0] = y0;
        for (int i = 0; i < Y.Length - 1; i++)
        {
            X[i] = a + i * h;
            var k1 = h * f(X[i], Y[i]);
            var k2 = h * f(X[i] + h / 2, Y[i] + k1 / 2);
            var k3 = h * f(X[i] + h / 2, Y[i] + k2 / 2);
            var k4 = h * f(X[i] + h, Y[i] + k3);
            Y[i + 1] = Y[i] + (k1+2*k2+2*k3+k4)/6;


            if (DisplayMode)
            {
                //Console.Write($"k1 = {k1} ");

                //Console.Write($"k2 = {k2} ");

                //Console.Write($"k3 = {k3} ");

                //Console.WriteLine($"k4 = {k4} ");

                Console.Write("(" + X[i] + ";" + Y[i] + ") ");
            }
        }
        RK[0] = X;
        RK[1] = Y;
        return RK;
    }

    static double[][] Runge_Kutta_2nd_Order(function3 f1, function3 f2, double a, double b, int N, double y0, double y1, bool DisplayMode = false)
    {
        double[][] RK = new double[3][];
        double h = (b - a) / N;
        double[] X = new double[N + 1];
        double[] Y0 = new double[N + 1];
        double[] Y1 = new double[N + 1];
        Y0[0] = y0;
        Y1[0] = y1;
        for (int i = 0; i < Y1.Length - 1; i++)
        {
            X[i] = a + i * h;
            var k1 = h * f1(X[i]        , Y0[i]          , Y1[i]);
            var m1 = h * f2(X[i]        , Y0[i]          , Y1[i]);
            var k2 = h * f1(X[i] + h / 2, Y0[i] + k1 / 2 , Y1[i] + m1 / 2);
            var m2 = h * f2(X[i] + h / 2, Y0[i] + k1 / 2 , Y1[i] + m1 / 2);
            var k3 = h * f1(X[i] + h / 2, Y0[i] + k2 / 2 , Y1[i] + m2 / 2);
            var m3 = h * f2(X[i] + h / 2, Y0[i] + k2 / 2 , Y1[i] + m2 / 2);
            var k4 = h * f1(X[i] + h    , Y0[i] + k3     , Y1[i] + m3);
            var m4 = h * f2(X[i] + h    , Y0[i] + k3     , Y1[i] + m3);
            
            Y0[i + 1] = Y0[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            Y1[i + 1] = Y1[i] + (m1 + 2 * m2 + 2 * m3 + m4) / 6;


            if (DisplayMode)
            {
                //Console.Write($"k1 = {k1} ");

                //Console.Write($"k2 = {k2} ");

                //Console.Write($"k3 = {k3} ");

                //Console.WriteLine($"k4 = {k4} ");

                Console.Write("("+X[i] + ";" + Y0[i]+") ");
            }
        }
        RK[0] = X;
        RK[1] = Y0;
        RK[2] = Y1;
        return RK;
    }






    static void DrawPolynom(double[] a)
    {
        for (int i = 0; i < a.Length-1; i++)
        {
            Console.Write($"{a[i]}*x^{i} + ");
        }
        Console.WriteLine($"{a[a.Length-1]}*x^{a.Length-1}");
    }
    //Интегралы
    static public double Integral(function f, double a, double b, int n, methodName methodName)
    {
        double h = (b - a) / n;
        double x = a;
        double S = 0.0;
        for (int i = 0; i < n; i++)
        {
            switch (methodName)
            {
                case (methodName.Left): { S += f(x); break; }
                case (methodName.Right): { S += f(x + h); break; }
                case (methodName.Middle): { S += f(x + h / 2); break; }
            }

            x += h;
        }
        return S * h;
    } //Метод Правых/Левых/Средних Прямоугольников
    public enum methodName
    {
        Left,
        Right,
        Middle
    } // Перечисление для метода Правых/Левых/Средних Прямоугольников
    static public double TrapezoidIntegral(function f, double a, double b, int n)
    {
        var S = 0.0;
        var h = (b - a) / n;
        var x = a;
        for (int i = 0; i < n; i++)
        {
            x = x + h;
            S += f(x) + f(x + h);
        }
        S = h / 2 * (S + f(a));
        return S;
    } //Метод трапеций

    static public double SimpsonIntegral(function f, double a, double b, int n)
    {
        if (n / 2 != 0) { n++; }
        var S = 0.0;
        var h = (b - a) / n;
        var c = 1;
        var x = a;
        for (int i = 0; i < n; i++)
        {
            x = x + h;
            S += (c + 3) * f(x);
            c = -c;
        }
        S = h / 3 * (S + f(a) + f(b));
        return S;
    } // Метод Симпсона




    //
    static public void DrawMatrix(double[,] matr)
    {
        int Rows = matr.GetLength(0);
        int Column = matr.GetLength(1);
        for (int i = 0; i < Rows; i++)
        {
            for (int j = 0; j < Column; j++)
            {
                Console.Write(String.Format("{0:0.000}", matr[i, j]) + " ");
            }
            Console.WriteLine();
        }
    } // Вспомогательные функции для отображения матриц
    static public void DrawMatrix(double[] matr)
    {
        int N = matr.Length;
        for (int i = 0; i < N; i++)
        {
            Console.Write(String.Format("{0:0.000}", matr[i]) + " ");
        }
        Console.WriteLine();
    }

    public static double[] CopyArray(double[] arr)
    {
        double[] copyied = new double[arr.Length];
        for (int i = 0; i < arr.Length; i++)
        {
            copyied[i] = arr[i];
        }
        return copyied;
    } // Удобное копирование массива

    //Решение СЛАУ
    static public void swapRows(ref double[,] matr, int row1, int row2)
    {
        var n = matr.GetLength(1);
        for (int i = 0; i < n; i++)
        {
            var temp = matr[row1, i];
            matr[row1, i] = matr[row2, i];
            matr[row2, i] = temp;
        }

    }

    static public double[] Gauss(double[,] matr)
    {

        int Rows = matr.GetLength(0);
        int Column = matr.GetLength(1);
        int i, k, q;
        double v;
        double[] answer = new double[Rows];

        for (q = 0; q < Rows; q++)
        {
            for (k = q + 1; k < Rows; k++)
            { //Находим самое большое по модулю число на диогонале для текущей строки 
                if (Math.Abs(matr[k, q]) > Math.Abs(matr[q, q]))
                {
                    swapRows(ref matr, k, q);
                }
            }
            //делаем главную диогональ единицами
            v = matr[q, q];//для этого делим на коэффицент диоганального элемента
            for (k = 0; k < Column; k++)
                matr[q, k] /= v;
            //обнуляем числа под единицами главной диогoнали
            for (i = q + 1; i < Rows; i++)
            {
                v = matr[i, q];//Для этого вычитаем из нижней строчки текущюю умноженуюю на значение коэффицента ненулевого элемента под диоганалью
                for (k = q; k < Column; k++)
                    matr[i, k] = matr[i, k] - matr[q, k] * v;
            }
        }// Получили диоганальную матрицу
        for (q = 0; q < Rows; q++)
            for (i = 0; i < (Rows - 1) - q; i++)
            {
                v = matr[i, (Column - 1) - q - 1];//элемент текущей строки соответствующий текущему свободному элементу (идем с конца тк последний элемент известен)
                for (k = Column - 1 - q - 1; k < Column; k++)
                    matr[i, k] = matr[i, k] - matr[(Rows - 1) - q, k] * v;//вычитаем из строки строку со свободным элементом умноженную на V
            }//при этом последнюю строку мы не трогаем и автоматически в ней будет результат
        for (i = 0; i < Rows; i++)
            answer[i] = matr[i, Column - 1];
        return answer;
    } // Метод Гаусса
    static public double[] iterativeMethod(double[,] matr, double[] b, double eps)
    {
        int counter = 0;
        (matr, b) = GetRightForm(matr, b);
        double[] X1;
        double[] X2 = CopyArray(b);
        do
        {

            X1 = CopyArray(X2);
            for (int i = 0; i <= matr.GetUpperBound(0); i++)
            {
                X2[i] = b[i];
                for (int j = 0; j <= matr.GetUpperBound(1); j++)
                {
                    X2[i] += matr[i, j] * X1[j];
                }
            }
            counter++;
            DrawMatrix(X2);
        }
        while (LeaveCondition(X1, X2) > eps);

        Console.WriteLine("Количество итераций = " + counter);
        return X2;
    }

    public static double LeaveCondition(double[] X1, double[] X2)
    {
        int N = X1.Length;
        double[] Temp = new double[N];
        for (int i = 0; i < N; ++i)
        {
            Temp[i] = Math.Abs(X1[i] - X2[i]);
        }
        double Max = Temp[0];
        for (int i = 1; i < N; ++i)
        {
            if (Max < Temp[i])
            {
                Max = Temp[i];
            }
        }
        return Max;
    }
    public static (double[,], double[]) GetRightForm(double[,] A, double[] b)
    {
        int N = b.Length;
        double[,] B = new double[N, N];
        double[] e = new double[N];
        for (int i = 0; i < N; ++i)
        {
            e[i] = b[i] / A[i, i];
            B[i, i] = 0;
        }
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                if (i == j)
                {
                    continue;
                }
                else
                {
                    B[i, j] = -(A[i, j] / A[i, i]);
                }
            }
        }
        return (B, e);
    }

    static public double[] zeidelMethod(double[,] matr, double[] b, double eps)
    {
        int counter = 0;
        (matr, b) = GetRightForm(matr, b);
        double[] X1 = CopyArray(b);
        double[] X2 = new double[X1.Length];
        do
        {
            X1 = CopyArray(X2);
            for (int i = 0; i <= matr.GetUpperBound(0); i++)
            {
                X2[i] = b[i];
                for (int j = 0; j <= matr.GetUpperBound(1); j++)
                {
                    if (j < i)
                    {
                        X2[i] += matr[i, j] * X2[j];
                    }
                    else
                    {
                        X2[i] += matr[i, j] * X1[j];
                    }
                }

            }
            counter++;

            DrawMatrix(X1);
        }
        while (LeaveCondition(X1, X2) > eps);
        return X2;
    }



    //Интерполяция
    public static double[] GetInterpolationCoefficient(double[] vX, double[] vY)
    {
        if (vX.Length != vY.Length)
        {
            Console.WriteLine("Длина входных данных: " + vX.Length);
            Console.WriteLine("Длина выходных данных: " + vY.Length);
            Console.WriteLine("...");
            return null;
        }
        double[,] matr = new double[vX.Length, vX.Length + 1];
        for (int i = 0; i < vX.Length; i++)
        {
            for (int j = 0; j < vX.Length; j++)
            {
                matr[i, j] = Math.Pow(vX[i], j);
            }
        }
        for (int i = 0; i < vX.Length; i++)
        {
            matr[i, matr.GetUpperBound(1)] = vY[i];
        }
        return Gauss(matr);
    }
    public static double GetValueFromPolynomial(double[] a, double x)
    {
        double y = 0;
        for (int i = 0; i < a.Length; i++)
        {
            //Console.Write(a[i]+"*" +"x^"+i+" ");
            y += a[i] * Math.Pow(x, i);
        }
        return y;
    }
    //Лагранж
    public static double LagrangeInterpolation(double[] vX, double[] vY, double x, bool DisplayMode = false)
    {
        double y = 0;
        for (int i = 0; i < vX.Length; i++)
        {
            double temp = vY[i];
            for (int j = 0; j < vX.Length; j++)
            {
                if (i == j) continue;
                temp *= x - vX[j];
                temp /= vX[i] - vX[j];
            }
            if (DisplayMode) Console.WriteLine("L" + i + " = " + temp);
            y += temp;
        }
        return y;
    }
    //Ньютон
    public static double NewtonInterpolation(double[] vX, double[] vY, double x, bool DisplayMode = false)
    {
        if (vX.Length != vY.Length)
        {
            Console.WriteLine("Длина входных данных: " + vX.Length);
            Console.WriteLine("Длина выходных данных: " + vY.Length);
            Console.WriteLine("Не совпадает");
            return 0;
        }

        int n = vX.Length;
        double h = vX[1] - vX[0];
        double[,] matr = new double[n, n];
        for (int i = 0; i < n; i++)
        {
            matr[i, 0] = vY[i];
        }
        for (int i = 1; i < n; i++)
        {
            for (int j = 0; j < n - i; j++)
            {
                matr[j, i] = matr[j + 1, i - 1] - matr[j, i - 1];
            }
        }



        if (DisplayMode)
        {
            Console.WriteLine("Таблица конечных разностей: ");
            DrawMatrix(matr);
        }
        double y = 0;

        double t = (x - vX[0]) / h;
        for (int i = 0; i < n; i++)
        {

            y += matr[0, i] * factorial(t, i) / factorial(i);
        }
        return y;
    }
    public static double NewtonInterpolationReversed(double[] vX, double[] vY, double y, bool DisplayMode = false)
    {
        if (vX.Length != vY.Length)
        {
            Console.WriteLine("Длина входных данных: " + vX.Length);
            Console.WriteLine("Длина выходных данных: " + vY.Length);
            Console.WriteLine("Не совпадает");
            return 0;
        }

        int n = vX.Length;
        double h = vX[1] - vX[0];

        double[,] matr = new double[n, n];
        for (int i = 0; i < n; i++)
        {
            matr[i, 0] = vY[i];
        }
        for (int i = 1; i < n; i++)
        {
            for (int j = 0; j < n - i; j++)
            {
                matr[j, i] = matr[j + 1, i - 1] - matr[j, i - 1];
            }
        }

        if (DisplayMode)
        {
            Console.WriteLine("Таблица конечных разностей: ");
            DrawMatrix(matr);
        }
        double x = 0.0;
        double x_old = 0;
        double t_ = (y - vY[0]) / (matr[0, 1]);
        do
        {
            x_old = x;
            double t = matr[0, 0] - y;
            for (int i = 2; i < n; i++)
            {
                t += matr[0, i] * factorial(t_, i) / factorial(i);
            }
            t /= -matr[0, 1];
            t_ = t;
            x = h * t_ + vX[0];
        } while (Math.Abs(x_old - x) > 0.000001);
        return x;
    }//Обратное
    public static double factorial(double n, int times = -1)
    {
        if (times == -1) times = (int)n;
        double temp = 1;
        while (times > 0)
        {
            temp *= n;
            n--;
            times--;
        }
        return temp;
    }


    //Сплайн
    public static double diff(function f, double x)
    {
        var h = 0.0001;
        double y_0 = (-3 * f(x - h) + 4 * f(x) - f(x + h)) / (2 * h);
        double y_1 = (-f(x - h) + f(x + h)) / (2 * h);
        double y_2 = (f(x - h) - 4 * f(x) + 3 * f(x + h)) / (2 * h);
        double y_ = (y_0 + y_1 + y_2) / 3;
        return y_;
    }//производная 1 порядка
    public static double diff_2(function f, double x)
    {
        var h = 0.0001;
        double y_0_ = diff(f, x - h);
        double y_1_ = diff(f, x);
        double y_2_ = diff(f, x + h);

        double y_0 = (-3 * y_0_ + 4 * y_1_ - y_2_) / (2 * h);
        double y_1 = (-y_0_ + y_2_) / (2 * h);
        double y_2 = (y_0_ - 4 * y_1_ + 3 * y_2_) / (2 * h);
        double y_ = (y_0 + y_1 + y_2) / 3;
        return y_;
    }//производная 2 порядка      //чето я забыл зачем они мне

    public static double spline(double[] vx, double[] vy, double x, bool DisplayMode = false)
    {

        var coeff = splineCoeff(vx, vy, DisplayMode);
        int n = vx.Length;
        double y = 0;
        for (int i = 0; i < n - 1; i++)
        {
            if (vx[i] <= x && vx[i + 1] >= x)
            {
                y = coeff[4 * i] + coeff[4 * i + 1] * (x - vx[i]) + coeff[4 * i + 2] * (x - vx[i]) * (x - vx[i]) + coeff[4 * i + 3] * (x - vx[i]) * (x - vx[i]) * (x - vx[i]);
                return y;
            }
        }
        Console.WriteLine("Выберите точку внутри узлов");
        return 0;
    }
    public static double[] splineCoeff(double[] vx, double[] vy, bool DisplayMode = false)
    {
        int n = vx.Length - 1;
        double[,] slau = new double[4 * n, 4 * n + 1];
        for (int i = 0; i < 4 * n; i += 1)
        {
            for (int j = 0; j <= 4 * n; j += 1)
            {
                var k = i % n;
                var h = vx[k + 1] - vx[k];
                if (i < n)
                {
                    if (j == 4 * n)
                    {
                        slau[i, j] = vy[i];
                    }
                    else if (j == 4 * i)
                    {
                        slau[i, j] = 1;
                    }
                    else
                    {
                        slau[i, j] = 0;
                    }
                }
                if (n <= i && i < 2 * n)
                {
                    if (j == 4 * k)
                    {
                        slau[i, j] = 1;
                    }
                    else if (j == 4 * k + 1)
                    {
                        slau[i, j] = h;
                    }
                    else if (j == 4 * k + 2)
                    {
                        slau[i, j] = h * h;
                    }
                    else if (j == 4 * k + 3)
                    {
                        slau[i, j] = h * h * h;
                    }
                    else if (j == 4 * n)
                    {
                        slau[i, j] = vy[k + 1];
                    }
                    else
                    {
                        slau[i, j] = 0;
                    }
                }
                if (2 * n <= i && i < 3 * n - 1)
                {
                    if (j == 4 * k + 1)//b
                    {
                        slau[i, j] = 1;
                    }
                    else if (j == 4 * k + 2)//c
                    {
                        slau[i, j] = 2 * h;
                    }
                    else if (j == 4 * k + 3)//d
                    {
                        slau[i, j] = 3 * h * h;
                    }
                    else if (j == 4 * k + 5)//b +1
                    {
                        slau[i, j] = -1;
                    }
                    else if (j == 4 * n)//s(x)
                    {
                        slau[i, j] = 0;
                    }
                    else
                    {
                        slau[i, j] = 0;
                    }
                }
                if (i == 3 * n - 1)
                {
                    if (j == 2)
                    {
                        slau[i, j] = 1;
                    }
                    else
                    {
                        slau[i, j] = 0;
                    }
                }
                if (3 * n <= i && i < 4 * n - 1)
                {

                    if (j == 4 * k + 2)//c
                    {
                        slau[i, j] = 2;
                    }
                    else if (j == 4 * k + 3)//d
                    {
                        slau[i, j] = 6 * h;
                    }
                    else if (j == 4 * k + 6)//c +1
                    {
                        slau[i, j] = -2;
                    }
                    else if (j == 4 * n)//s(x)
                    {
                        slau[i, j] = 0;
                    }
                    else
                    {
                        slau[i, j] = 0;
                    }
                }
                if (i == 4 * n - 1)
                {
                    if (j == 4 * n - 1)//d n
                    {
                        slau[i, j] = 3 * h;
                    }
                    else if (j == 4 * n - 2)//d n
                    {
                        slau[i, j] = 1;
                    }
                    else
                    {
                        slau[i, j] = 0;
                    }
                }
            }
        }

        var coeff = Gauss(slau);

        if (DisplayMode)
        {
            Console.WriteLine("Система линейных уравнений: ");
            DrawMatrix(slau);
            Console.WriteLine("Коэффиценты: ");
            DrawMatrix(coeff);
        }
        Console.WriteLine("Выберите точку внутри узлов");
        return coeff;
    }//Получить коэфиценты сплайна

    //Апроксимация

    public static double SumSeq(double[] sequence)
    {
        double temp = 0;
        for (int i = 0; i < sequence.Length; i++)
        {
            temp += sequence[i];
        }
        return temp;
    }

    public static double[] SeqPow(double[] sequence, int Power)
    {
        double[] temp = new double[sequence.Length];
        for (int i = 0; i < temp.Length; i++)
        {
            temp[i] = Math.Pow(sequence[i], Power);
        }
        return temp;
    }

    public static double[] AproximationPolynom(double[] vx, double[] vy, int power, bool DisplayMode = false)
    {
        int n = vx.Length;
        if (n <= power)
        {
            Console.WriteLine("Степень должна быть меньше количества исходных данных");
        }
        double[,] matr = new double[power + 1, power + 2];
        for (int i = 0; i < power + 1; i++)
        {
            for (int j = 0; j < power + 1; j++)
            {
                matr[i, j] = SumSeq(SeqPow(vx, i + j));
            }
        }
        for (int j = 0; j < power + 1; j++)
        {
            double[] xy = new double[n];
            double[] x_pow = SeqPow(vx, j);
            for (int i = 0; i < vx.Length; i++)
            {
                xy[i] = x_pow[i] * vy[i];
            }
            
            matr[j, power + 1] = SumSeq(xy);
        }
        if (DisplayMode) DrawMatrix(matr);


        return Gauss(matr);

    }

    public static double[] AproximationExp(double[] vx, double[] vy, bool DisplayMode = false)
    {
        double[] new_vy = new double[vx.Length];
        for (int i = 0; i < new_vy.Length; i++)
        {
            new_vy[i] = Math.Log(vy[i]);
        }

        double[] temp = AproximationPolynom(vx, new_vy, 1,true);
        DrawMatrix(temp);
        double[] answer = new double[] { Math.Exp(temp[0]), temp[1] };
        if (DisplayMode) Console.WriteLine($"{answer[1]}*exp(x*{answer[0]})");
        return answer;

    }

    public static double[] AproximationLn(double[] vx, double[] vy, bool DisplayMode = false)
    {
        double[] new_vx = new double[vx.Length];
        for (int i = 0; i < new_vx.Length; i++)
        {
            new_vx[i] = Math.Log(vx[i]);
        }
        double[] answer = AproximationPolynom(new_vx, vy, 1);
        if (DisplayMode) Console.WriteLine($"{answer[1]}*Ln(x)+{answer[0]}");
        return answer;

    }

    public static double[] AproximationPow(double[] vx, double[] vy, bool DisplayMode = false)
    {
        double[] new_vx = new double[vx.Length];
        double[] new_vy = new double[vx.Length];
        for (int i = 0; i < new_vx.Length; i++)
        {
            new_vx[i] = Math.Log(vx[i]);
            new_vy[i] = Math.Log(vy[i]);
        }
        double[] temp = AproximationPolynom(new_vx, new_vy, 1);
        double[] answer = new double[] { Math.Exp(temp[0]), temp[1] };
        if (DisplayMode) Console.WriteLine($"{answer[0]}*x^{answer[1]}");
        return answer;

    }

    public static double[] AproximationGyperbolic(double[] vx, double[] vy, bool DisplayMode = false)
    {
        double[] new_vx = new double[vx.Length];
        for (int i = 0; i < new_vx.Length; i++)
        {
            new_vx[i] = 1 / (vx[i]);
        }
        double[] answer = AproximationPolynom(new_vx, vy, 1);
        if (DisplayMode) Console.WriteLine($"{answer[1]}*(1/x)+{answer[0]}");
        return answer;

    }

    public static double Corr(double[] vxx, double[] vyy, function function)
    {
        double[] vx = CopyArray(vxx);
        double[] vy = CopyArray(vyy);
        int N = vy.Length;
        double fx_sr = 0;
        double vy_sr = 0;
        double[] fx = new double[N];
        for (int i = 0; i < N; i++)
        {
            fx[i] = function(vx[i]);
            fx_sr += fx[i];
            vy_sr += vy[i];
        }
        fx_sr /= N;
        vy_sr /= N;
        double[] fx_vy = new double[N];
        for (int i = 0; i < N; i++)
        {
            fx[i] -= fx_sr;
            vy[i] -= vy_sr;
            fx_vy[i] = fx[i] * vy[i];
        }
        double r = (SumSeq(fx_vy)) / ((Math.Sqrt(SumSeq(SeqPow(fx, 2)))) * (Math.Sqrt(SumSeq(SeqPow(vy, 2)))));
        return r;
    }
}