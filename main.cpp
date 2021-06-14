#include <iostream>
#include <vector>
#include <iomanip>
#include "math.h"
#define GNUPLOT_NAME "gnuplot -persist"


using namespace std;
static vector<double> l;
static double epsilon = 0.0001;

vector<vector<double>> operator+(vector<vector<double>> matrixA, vector<vector<double>> matrixB);

vector<double> operator+(vector<double> matrixA, vector<double> matrixB);

vector<vector<double>> operator-(vector<vector<double>> matrixA, vector<vector<double>> matrixB);

vector<vector<double>> operator*(vector<vector<double>> matrixA, vector<vector<double>> matrixB);

vector<double> operator*(vector<vector<double>> matrixA, vector<double> matrixB);

vector<double> operator*(vector<double> matrixA, double coeff);

bool operator!=(vector<vector<double>> matrixA, vector<vector<double>> matrixB);

double findMax(vector<double> a) {
    int max = 0;
    int k = 0;
    for (int i = 0; i < a.size(); ++i) {
        if (a[i] * a[i] > max * max) {
            max = a[i];
            k = i;
        }
    }
    return k;
}

class Matrix {
public:
    void setN(int n) {
        Matrix::n = n;
    }

    void setM(int m) {
        Matrix::m = m;
    }

    static vector<double> doMultByC(vector<double> line, double coeff) {
        for (int i = 0; i < line.size(); ++i) {
            line[i] = line[i] * coeff;
        }
        return line;
    }

    static vector<double> doAddLines(vector<double> l1, vector<double> l2) {
        vector<double> res;
        for (int i = 0; i < l1.size(); ++i) {
            res.push_back(l1[i] + l2[i]);
        }
        return res;
    }

    /**
     * addition of two matrices
     * checks whether they have the same amount of rows and columns
     * @param matrixA first matrix
     * @param matrixB second matrix
     * @return resulting matrix
     */
    static vector<vector<double>> doAddition(vector<vector<double>> matrixA, vector<vector<double>> matrixB) {
        vector<vector<double>> matrixC;

        if (matrixA.size() == matrixB.size() && matrixA[0].size() == matrixB[0].size()) {
            for (int i = 0; i < matrixA.size(); ++i) {
                vector<double> row;
                for (int j = 0; j < matrixA[0].size(); ++j) {
                    row.push_back(matrixA[i][j] + matrixB[i][j]);
                }
                matrixC.push_back(row);
            }
        } else {
            cout << "Error: the dimensional problem occurred" << endl;
        }
        return matrixC;
    }

    /**
     * subtracting two matrices
     * checks whether they have the same amount of rows and columns
     * @param matrixA first matrix
     * @param matrixB second matrix
     * @return resulting matrix
     */
    static vector<vector<double>> doSubtruction(vector<vector<double>> matrixA, vector<vector<double>> matrixB) {
        vector<vector<double>> matrixC;

        if (matrixA.size() == matrixB.size() && matrixA[0].size() == matrixB[0].size()) {
            for (int i = 0; i < matrixA.size(); ++i) {
                vector<double> row;
                for (int j = 0; j < matrixA[0].size(); ++j) {
                    row.push_back(matrixA[i][j] - matrixB[i][j]);
                }
                matrixC.push_back(row);
            }
        } else {
            cout << "Error: the dimensional problem occurred" << endl;
        }
        return matrixC;
    }


    /**
     * multiplying two matrices
     * checks whether amount of columns of the first one is equal to amount of rows of the second one
     * @param matrixA
     * @param matrixB
     * @return resulting matrix
     */
    static vector<vector<double>> doMultiplication(vector<vector<double>> matrixA, vector<vector<double>> matrixB) {
        vector<vector<double>> matrixC;
        if (matrixA[0].size() == matrixB.size()) {
            for (int i = 0; i < matrixA.size(); ++i) {
                double res = 0;
                int k = 0;
                vector<double> row;
                while (k != matrixB[0].size()) {
                    for (int j = 0; j < matrixA[0].size(); ++j) {
                        res += (matrixA[i][j] * matrixB[j][k]);
                    }
                    k++;
                    row.push_back(res);
                    res = 0;
                }
                matrixC.push_back(row);
            }
        } else {
            cout << "Error: the dimensional problem occurred" << endl;
        }

        return matrixC;
    }


    /**
     * transposing given matrix
     * @param matrixA
     * @return
     */
    static vector<vector<double>> doTransposeMatrix(vector<vector<double>> matrixA) {
        vector<vector<double>> matrixT;
        for (int i = 0; i < matrixA[0].size(); ++i) {
            vector<double> row;
            for (int j = 0; j < matrixA.size(); ++j) {
                row.push_back(matrixA[j][i]);
            }
            matrixT.push_back(row);
            row.clear();
        }
        return matrixT;
    };

    static bool isEqual(vector<vector<double>> matrixA, vector<vector<double>> matrixB) {
        vector<vector<double>> A = matrixA;
        vector<vector<double>> B = matrixB;
        if (A.size() == B.size()) {
            for (int i = 0; i < A.size(); ++i) {
                for (int j = 0; j < A.size(); ++j) {
                    if (A[i][j] != B[i][j]) {
                        return true;
                    }
                }
            }
            return false;
        } else {
            return true;
        }
    }


    static vector<vector<double>> insertingMatrix() {
        int n, m;
        cin >> n >> m;
        vector<vector<double>> matrixA;
        for (int i = 0; i < n; ++i) {
            vector<double> row;
            for (int j = 0; j < m; ++j) {
                int element;
                cin >> element;
                row.push_back(element);
            }
            matrixA.push_back(row);
        }
        return matrixA;
    }

    static void displayMatrix(vector<vector<double>> res) {
        for (int i = 0; i < res.size(); ++i) {
            for (int j = 0; j < res[0].size(); ++j) {
                if (abs(res[i][j]) <= epsilon) {
                    cout << fixed << setprecision(4) << (double) 0 << " ";
                } else {
                    cout << fixed << setprecision(4) << res[i][j] << " ";
                }
            }
            cout << endl;
        }
    }


private:
    int n;//number of rows
    int m;//number of columns
};

class PermutationMatrix;

class ColumnVector;

class IdentityMatrix;

class SquareMatrix : public Matrix {
public:

    static vector<vector<double>> insertingMatrix() {
        int n;
        cin >> n;
        vector<vector<double>> matrixA;
        for (int i = 0; i < n; ++i) {
            vector<double> row;
            for (int j = 0; j < n; ++j) {
                double element;
                cin >> element;
                row.push_back(element);
            }
            matrixA.push_back(row);
        }
        return matrixA;
    }


};

class IdentityMatrix : public SquareMatrix {
public:

    static vector<vector<double>> createIdentity(int n) {
        vector<vector<double>> identity;
        if (n != 0) {
            for (int i = 0; i < n; ++i) {
                vector<double> row;
                for (int j = 0; j < n; ++j) {
                    if (i == j) {
                        row.push_back(1);
                    } else {
                        row.push_back(0);
                    }
                }
                identity.push_back(row);
            }
        }
        return identity;
    }

    static int n;
};

class EliminationMatrix;

static vector<EliminationMatrix> elimMatrixArray;

class PermutationMatrix : public IdentityMatrix {
public:
    static vector<vector<double>> createPermutation(int n, double r1, double r2) {
        vector<vector<double>> identity;
        if (n != 0) {
            for (int i = 0; i < n; ++i) {
                vector<double> row;
                for (int j = 0; j < n; ++j) {
                    if (i == j) {
                        row.push_back(1);
                    } else {
                        row.push_back(0);
                    }
                }
                identity.push_back(row);
            }
        }

        vector<double> temp = identity[r2 - 1];
        identity[r2 - 1] = identity[r1 - 1];
        identity[r1 - 1] = temp;

        return identity;
    }

    static vector<vector<double>> doPermutation(vector<vector<double>> ma, int c) {
        vector<vector<double>> matrix = ma;
        vector<double> col = Matrix::doTransposeMatrix(matrix)[c - 1];
        int row = findMax(col);
        return createPermutation(matrix.size(), row + 1, 1) * matrix;
    }
};

class EliminationMatrix : public SquareMatrix {
public:
    int r;
    int c;
    vector<vector<int>> elimMatrix;

    static vector<vector<double>> swap_row(vector<vector<double>> mat, int i, int j) {
        vector<double> temp = mat[i];
        mat[i] = mat[j];
        mat[j] = temp;
        return mat;
    }

    EliminationMatrix(int r, int c, const vector<vector<int>> &elimMatrix) : r(r), c(c), elimMatrix(elimMatrix) {}

    static void insert(int coefficient, int r, int c, int n) {
        vector<vector<int>> elim;
        for (int i = 0; i < n; ++i) {
            vector<int> row;
            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    row.push_back(1);
                } else if (i == r - 1 && j == c - 1) {
                    row.push_back(coefficient);
                } else {
                    row.push_back(0);
                }
            }
            elim.push_back(row);
        }
        EliminationMatrix *instance = new EliminationMatrix(r, c, elim);
        elimMatrixArray.push_back(*instance);
    }

    static vector<vector<double>> doUmatrix(vector<vector<double>> vec) {
        int counter = 1;
        int n = vec.size();
        vector<vector<double>> matrix = vec;
        for (int k = 0; k < n; k++) {
            int i_max = k;
            int v_max = matrix[i_max][k];

            for (int i = k + 1; i < n; i++)
                if (abs(matrix[i][k]) > abs(v_max))
                    v_max = matrix[i][k], i_max = i;

            if (i_max != k) {
                matrix = swap_row(matrix, k, i_max);
                cout << "step #" << counter << ": permutation" << endl;
                Matrix::displayMatrix(matrix);
            }


            for (int i = k + 1; i < n; i++) {
                double f = matrix[i][k] / matrix[k][k];

                for (int j = k + 1; j <= n; j++) {
                    matrix[i][j] -= matrix[k][j] * f;
                    EliminationMatrix::insert(-f, counter + 1, 1, 3);
                }

                matrix[i][k] = 0;
                cout << "step #" << counter << ": elimination" << endl;
                Matrix::displayMatrix(matrix);
                counter++;
            }
        }
        return matrix;
    }

    static double calcDeterminant(vector<vector<double>> mat) {
        double det = 1;
        for (int i = 0; i < mat.size(); ++i) {
            for (int j = 0; j < mat.size(); ++j) {
                if (i == j) {
                    det = det * mat[i][j];
                }
            }
        }
        return det;
    }

};

class InverseMatrix : public EliminationMatrix {
public:
    static vector<vector<double>> doAugemented(vector<vector<double>> matrix) {
        vector<vector<double>> newMat;
        for (int i = 0; i < matrix.size(); ++i) {
            vector<double> row;
            for (int j = 0; j < 2 * matrix.size(); ++j) {
                if (j < matrix.size()) {
                    row.push_back(matrix[i][j]);
                } else if (j == matrix.size() + i) {
                    row.push_back(1);
                } else {
                    row.push_back(0);
                }
            }
            newMat.push_back(row);
        }
        return newMat;
    }

    static vector<vector<double>> doInverse(vector<vector<double>> matrix) {
        int counter = 1;
        int rem = counter;
        int n = matrix.size();
        for (int k = 0; k < n; k++) {
            int i_max = k;
            int v_max = matrix[i_max][k];

            rem = counter;

            for (int h = k; h < n; h++) {
                for (int i = h + 1; i < n; i++) {
                    if (abs(matrix[i][h]) > abs(v_max)) {
                        v_max = matrix[i][h], i_max = i;
                    }
                }
                if (i_max != h && i_max != k) {
                    matrix = swap_row(matrix, h, i_max);
                    //cout << "step #" << counter << ": permutation" << endl;
                    counter++;
//                    Matrix::displayMatrix(matrix);
                }
                if (rem == counter || matrix[k + 1][h] != 0) {
                    break;
                }
            }


            for (int i = k + 1; i < n; i++) {
                double f = matrix[i][k] / matrix[k][k];
                if (f != 0) {
                    // for (int j = k + 1; j <= n; j++) {
                    //    matrix[i][j] -= matrix[k][j] * f;
                    vector<double> row = matrix[k] * f * -1;
                    matrix[i] = matrix[i] + row;
                    //EliminationMatrix::insert(-1.00*f, counter + 1, 1, 3);
                    //  }
                    //  matrix[i][k] = 0;
//                    cout << "step #" << counter << ": elimination" << endl;
//                    Matrix::displayMatrix(matrix);
                    counter++;
                }
            }
        }

        double count = n - 1;
        //cout << "Way back:" << endl;
        for (int i = n - 1; i >= 1; --i) {
            for (int j = i; j > 0; --j) {
                double k = matrix[j - 1][i] / matrix[i][count];
                vector<double> row = matrix[i] * k * -1;
                matrix[j - 1] = matrix[j - 1] + row;
               // cout << "step #" << counter << ": elimination" << endl;
                counter++;
                //Matrix::displayMatrix(matrix);

            }
           count--;
        }
        return matrix;
    }

    static vector<vector<double>> normalize(vector<vector<double>> matrix) {
        vector<vector<double>> norm;
        int counter = 0;
        for (int i = 0; i < matrix.size(); ++i) {
            double coeff = 1 / matrix[i][counter];
            norm.push_back(matrix[i] * coeff);
            counter++;
        }
        return norm;
    }

    static vector<vector<double>> reversed(vector<vector<double>> matrix) {
        vector<vector<double>> norm;
        int counter = 0;
        for (int i = 0; i < matrix.size(); ++i) {
            vector<double> row;
            for (int j = 0; j < matrix.size(); ++j) {
                row.push_back(matrix[i][matrix.size() + j]);
            }
            norm.push_back(row);
        }
        return norm;
    }

};

class ColumnVector : InverseMatrix {
public:

    static vector<double> swap(vector<double> mat, int i, int j) {
        double temp = mat[i];
        mat[i] = mat[j];
        mat[j] = temp;
        return mat;
    }

    static void displayVec(vector<double> vec) {
        for (int i = 0; i < vec.size(); ++i) {
            if (abs(vec[i]) < epsilon) {
                cout << 0.00 << endl;
            } else {
                cout << vec[i] << endl;
            }
        }
    }

    static vector<double> insertingVector() {
        int n;
        cin >> n;
        vector<double> matrixA;
        for (int i = 0; i < n; ++i) {
            double element;
            cin >> element;
            matrixA.push_back(element);
        }
        return matrixA;
    }

    static vector<vector<double>> doUmatrix(vector<vector<double>> vec, vector<double> v) {
        int counter = 1;
        int n = vec.size();
        vector<vector<double>> matrix = vec;
        for (int k = 0; k < n; k++) {
            int i_max = k;
            double v_max = matrix[i_max][k];

            for (int i = k + 1; i < n; i++)
                if (abs(matrix[i][k]) - abs(v_max) > epsilon)
                    v_max = matrix[i][k], i_max = i;

            if (i_max != k) {
                matrix = swap_row(matrix, k, i_max);
                v = swap(v, k, i_max);
                cout << "step #" << counter << ": permutation" << endl;
                counter++;
                Matrix::displayMatrix(matrix);
                displayVec(v);
            }


            for (int i = k + 1; i < n; i++) {
                double f = matrix[i][k] / matrix[k][k];
                if (f != 0) {
                    vector<double> row = matrix[k] * f * -1;
                    matrix[i] = matrix[i] + row;
                    v[i] = v[i] + v[k] * f * -1;
                    cout << "step #" << counter << ": elimination" << endl;
                    Matrix::displayMatrix(matrix);
                    displayVec(v);
                    counter++;
                }
            }
        }
        double count = n - 1;
        cout << "Way back:" << endl;
        for (int i = n - 1; i >= 1; --i) {
            for (int j = i; j > 0; --j) {
                double k = matrix[j - 1][i] / matrix[i][count];
                vector<double> row = matrix[i] * k * -1;
                matrix[j - 1] = matrix[j - 1] + row;
                v[j - 1] = v[j - 1] + v[i] * k * -1;
                cout << "step #" << counter << ": elimination" << endl;
                counter++;
                Matrix::displayMatrix(matrix);
                displayVec(v);
            }
            count--;
        }
        l = v;
        return matrix;
    }

    static void normalize(vector<vector<double>> matrix, vector<double> v) {
        vector<vector<double>> norm;
        int counter = 0;
        for (int i = 0; i < matrix.size(); ++i) {
            double coeff = 1 / matrix[i][counter];
            l[i] = l[i] * coeff;
            norm.push_back(matrix[i] * coeff);
            counter++;
        }
        Matrix::displayMatrix(norm);
        displayVec(l);
    }

    static void displayMatrix(vector<vector<double>> res, vector<double> vec) {
        for (int i = 0; i < res.size(); ++i) {
            for (int j = 0; j < res[0].size(); ++j) {
                cout << fixed << setprecision(2) << res[i][j] << " ";
            }
            cout << endl;
        }
        displayVec(vec);
    }

};


vector<vector<double>> operator+(vector<vector<double>> matrixA, vector<vector<double>> matrixB) {
    return Matrix::doAddition(matrixA, matrixB);
}

vector<double> operator+(vector<double> matrixA, vector<double> matrixB) {
    return Matrix::doAddLines(matrixA, matrixB);
}

vector<vector<double>> operator-(vector<vector<double>> matrixA, vector<vector<double>> matrixB) {
    return Matrix::doSubtruction(matrixA, matrixB);
}

vector<vector<double>> operator*(vector<vector<double>> matrixA, vector<vector<double>> matrixB) {
    return Matrix::doMultiplication(matrixA, matrixB);
}

vector<double> operator*(vector<double> matrixA, double coeff) {
    return Matrix::doMultByC(matrixA, coeff);
}

vector<double> operator*(vector<vector<double>> matrixA, vector<double> matrixB) {
    vector<double> matrixC;
    if (matrixA[0].size() == matrixB.size()) {
        for (int i = 0; i < matrixA.size(); ++i) {
            double res = 0;
            for (int j = 0; j < matrixA[0].size(); ++j) {
                res += matrixA[i][j] * matrixB[j];
                //matrixC.push_back(matrixA[i][j]*matrixB[j]);
            }
            matrixC.push_back(res);
        }
    } else {
        cout << "Error: the dimensional problem occurred" << endl;
    }

    return matrixC;
}

double operator*(vector<double> a, vector<double> b){
    double res=0;
    for (int i = 0; i < a.size(); ++i) {
        res+= a[i]*b[i];
    }
    return res;
}

bool operator!=(vector<vector<double>> matrixA, vector<vector<double>> matrixB) {
    return Matrix::isEqual(matrixA, matrixB);
}

static void displayVec(vector<double> vec) {
    for (int i = 0; i < vec.size(); ++i) {
        if (abs(vec[i]) < epsilon) {
            cout << 0.0000 << endl;
        } else {
            cout<< fixed << setprecision(4)<< vec[i] << endl;
        }
    }
}

static vector<vector<double>> doBmatrix(vector<vector<double>> A) {
    vector<vector<double>> B;
    for (int i = 0; i < A.size(); ++i) {
        vector<double> row;
        for (int j = 0; j < A.size(); ++j) {
            if (i == j) {
                row.push_back(0);
            } else {
                row.push_back(-A[i][j] / A[i][i]);
            }
        }
        B.push_back(row);
    }
    return B;
}

static vector<double> doCvector(vector<vector<double>> A, vector<double> b) {
    vector<vector<double>> B = doBmatrix(A);
    vector<double> C;
    for (int i = 0; i < A.size(); ++i) {
        if (A[i][i] != 0) {
            C.push_back(b[i] / A[i][i]);
        } else {
            cout << "The method is not applicable!";
        }
    }
    return C;
}

static double delta(vector<double> x1, vector<double> x2) {
    double eps = 0;
    for (int i = 0; i < x1.size(); ++i) {
        eps += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }
    return sqrt(eps);
}

static bool checkJ(vector<vector<double>> matrix, vector<double> b, vector<double> x, double eps) {
    vector<double> C = doCvector(matrix, b);
    vector<vector<double>> B = doBmatrix(matrix);
    vector<double> x1 = C;
    vector<double> x2 = (B * x1) + C;
    double d = delta(x1, x2);
    double f = d;
    x1 = x2;
    x2 = (B * x2) + C;
    d = delta(x1, x2);
    if(d>f){
        return false;
    }else{
        return true;
    }
}

static void doJacobiMethod(vector<vector<double>> matrix, vector<double> b, vector<double> x, double eps) {
    if (checkJ(matrix, b, x, eps)) {
        int counter = 0;
        vector<double> C = doCvector(matrix, b);
        vector<vector<double>> B = doBmatrix(matrix);
        cout << "alpha: " << endl;
        Matrix::displayMatrix(B);
        cout << "beta:" << endl;
        displayVec(C);
        cout << "x(" << counter << "):" << endl;
        vector<double> x1 = C;
        displayVec(x1);
        vector<double> x2 = (B * x1) + C;
        double d = delta(x1, x2);
        while (d >= eps) {
            counter++;
            cout << "e: " << d << endl;
            cout << "x(" << counter << "):" << endl;
            displayVec(x2);
            x1 = x2;
            x2 = (B * x2) + C;
            d = delta(x1, x2);
        }
        counter++;
        cout << "e: " << d << endl;
        cout << "x(" << counter << "):" << endl;
        displayVec(x2);
    } else {
        cout << "The method is not applicable!"<<endl;
    }
}

static vector<vector<double>> doB(vector<vector<double>> matrix){
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < matrix.size(); ++j) {
            if(i<j){
                matrix[i][j]=0;
            }
        }
    }
    return matrix;
}
static vector<vector<double>> doC(vector<vector<double>> matrix){
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < matrix.size(); ++j) {
            if(i>j){
                matrix[i][j]=0;
            }
        }
    }
    return matrix;
}

static void doSeidelMethod(vector<vector<double>> matrix,vector<double> b, vector<double> x, double eps){
    if(checkJ(matrix, b, x,eps)) {
        int counter = 0;
        vector<double> beta = doCvector(matrix, b);
        vector<vector<double>> alpha = doBmatrix(matrix);
        cout << "beta:" << endl;
        displayVec(beta);
        cout << "alpha: " << endl;
        Matrix::displayMatrix(alpha);
        cout << "B:" << endl;
        vector<vector<double>> B = doB(alpha);
        Matrix::displayMatrix(B);
        cout << "C:" << endl;
        vector<vector<double>> C = doC(alpha);
        Matrix::displayMatrix(C);
        vector<vector<double>> identity = IdentityMatrix::createIdentity(B.size());
        cout << "I-B:" << endl;
        Matrix::displayMatrix(identity - B);
        vector<vector<double>> res = identity - B;
        res = InverseMatrix::doAugemented(res);
        res = InverseMatrix::doInverse(res);
        vector<vector<double>> rev;
        for (int i = 0; i < res.size(); ++i) {
            vector<double> row;
            for (int j = 0; j < res.size(); ++j) {
                row.push_back(res[i][res.size() + j]);
            }
            rev.push_back(row);
        }
        cout << "(I-B)_-1:" << endl;
        Matrix::displayMatrix(rev);
        cout << "x(" << counter << "):" << endl;
        displayVec(beta);
        vector<double> siVec = beta;
        double d = 10;
        while (d >= eps) {
            vector<double> was = beta;
            for (int i = 0; i < beta.size(); ++i) {
                beta[i] = alpha[i] * beta + siVec[i];
            }
            d = delta(beta, was);
            counter++;
            cout << "e: " << d << endl;
            cout << "x(" << counter << "):" << endl;
            displayVec(beta);

        }
    }else{
        cout<<"The method is not applicable!";
    }

}

static vector<vector<double>> insertLSq(){
    int n;
    cin>>n;
    vector<vector<double>> matrix;
    for (int i = 0; i < n; ++i) {
        vector<double> row;
        for (int j = 0; j < 2; ++j) {
            double el;
            cin>>el;
            row.push_back(el);
        }
        matrix.push_back(row);
    }
    return matrix;
}

static void doLSqpprox(vector<vector<double>> matrix, int deg){
    vector<vector<double>> A;
    vector<double> b;
    for (int i = 0; i < matrix.size(); ++i) {
        b.push_back(matrix[i][1]);
    }
    for (int i = 0; i < matrix.size(); ++i) {
        vector<double>row;
        for (int j = 0; j < deg+1; ++j) {
            if(j==0){
                row.push_back(1);
            }else{
                row.push_back(pow(matrix[i][0],j));
               // A[i][j]=pow(matrix[i][0],deg);
            }
        }
        A.push_back(row);
    }
    cout<<"A:"<<endl;
    //Matrix::displayMatrix(A);
    vector<vector<double>> A_T = A;
    A_T=SquareMatrix::doTransposeMatrix(A);
    cout<<"A_T*A:"<<endl;
    Matrix::displayMatrix(A_T*A);
    vector<vector<double>> res = A_T*A;
    res= InverseMatrix::doAugemented(res);
    //Matrix::displayMatrix(res);
    res = InverseMatrix::doInverse(res);
    //Matrix::displayMatrix(res);
    res = InverseMatrix::normalize(res);
    //Matrix::displayMatrix(res);
    vector<vector<double>> rev;
    for (int i = 0; i < res.size(); ++i) {
        vector<double> row;
        for (int j = 0; j < res.size(); ++j) {
            row.push_back(res[i][res.size() + j]);
        }
        rev.push_back(row);
    }
    cout<<"(A_T*A)^-1:"<<endl;
    Matrix::displayMatrix(rev);
    cout<<"A_T*b:"<<endl;
    displayVec(A_T*b);
    cout<<"x~:"<<endl;
    displayVec(rev*(A_T*b));
    vector<double> gnu = rev*(A_T*b);

    FILE *pipe = popen(GNUPLOT_NAME, "w");

    if (pipe != NULL) {
        const double pi = 3.14;
        const double npoints = 200;
        const double step = 4 * pi / npoints;

        fprintf(pipe, "%s\n", "plot '-' using 1:2 title 'exp' with points pointtype 5, '-' using 1:2 title 'appr' with lines");
        for (int i = 0; i < matrix.size(); ++i) {
            fprintf(pipe,  "%f\t%f\n", matrix[i][0], matrix[i][1]);
        }
        fprintf(pipe, "%s\n", "e");
        fprintf(pipe, "%s\n", "plot '-' using 1:2 title 'my graph' with lines");

        
        for (int i = 0; i < npoints + 1; ++i) {
            double x = -2 * pi + i * step+7;
            double y = 0;
            for (int j = gnu.size()-1; j >-1; --j) {
                y+=pow(x,j)*gnu[j];
            }
            fprintf(pipe, "%f\t%f\n", x, y);
        }
        fprintf(pipe,"\n");

        fprintf(pipe, "%s\n", "plot '-' using 1:2 title 'exp' with points pointtype 5, '-' using 1:2 title 'appr' with lines");


        fprintf(pipe, "%s\n", "e");
        fflush(pipe);

        pclose(pipe);
    } else {
        std::cout << "Could not open pipe!" << std::endl;
    }
}


int main() {
    vector<vector<double>> iden = insertLSq();
    int deg;
    cin>>deg;
    doLSqpprox(iden, deg);
    return 0;
}
