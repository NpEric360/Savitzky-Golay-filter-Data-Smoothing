// 4AX3_A2.cpp : This file contains the 'main' function. Program execution begins and ends there.


#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace::std;


/*
	Algorithm:
	1. Create a window of the original signal with size of m
	2. Build phi matrix given m and n
	3. Calculate a polynomial of degree n that fits the points in this window by building Theta with y values of that window
	4. Evaluate this polynomial at the midpoint of the window by multiply phi*theta;
		Note: due to centering; the first and last (m-1)/2 samples will be "lost" so for these edge cases:
		just take the first 3 for first sample, and last 3 for last sample
	5. Shift the window right by 1, calculate a new polynomial (theta by inputting the new corresponding y values).
	6. Repeat until end of file

*/

//Function 1: Compute the coefficients of the polynomial order n, m data points = sgolay
MatrixXd comp_coeff(MatrixXd phi, MatrixXd y, int n) {
	//calculate theta
	MatrixXd theta(n + 1, 1); //n+1 coefficients: i.e. linear requires 2 coefficients
	theta = (phi.transpose() * phi).inverse() * phi.transpose() * y;
	return theta; //theta
}

//Function 2: Apply the coefficients from comp_coeff as a filter to the data = sgolayfilt
//double* coeff_filter(double* coeff, int n, double* data, int m) {  
//	return 0; //return the filtered y, ****edge cases
//}

MatrixXd coeff_filter(int n, int m, MatrixXd phi, MatrixXd theta) {
	//if first 2 or last 2 items: take the first two values of phi*theta 
	MatrixXd y_filtered(m, 1);
	y_filtered = phi * theta;

	return y_filtered;
}


int main()
{
	//Part 0: Grab user input for polynomial order n, window size m:
	int m, n; //window size, polynomial order
	const int array_length = 2048;
	cout << "Input frame length m (must be odd!):";
	cin >> m;

	cout << "Input polynomial order n (must be less than m!):";
	cin >> n;
	//cout << "How many samples are in your dataset? Given = 2048: ____";
	//cin >> array_length;


//Part 1: Initialize Variables and matrices



	Eigen::MatrixXd y(m, 1);
	Eigen::MatrixXd phi(m, n + 1); //m rows, n+1 columns x^0 + ....x^n
	Eigen::MatrixXd theta(n + 1, 1);




	double data[array_length];


	//Part 2: Read .txt dataset
	std::ifstream file("fft_data.txt");
	std::string str;

	int i = 0;
	while (std::getline(file, str)) {
		//std::cout << str;
		data[i] = std::stod(str); //stod convert str to double
		cout << data[i] << endl;
		i++;
	}

	//Part 3: Write to .txt dataset FIX THIS
	ofstream myfile;
	myfile.open("a2_filtered_output_n2_m63.txt");
	//myfile << filtered_y;



//Part 4: Build Phi Vandermonde Matrix:
	//size  = m x (n+1)
	// centered about m
	//[x_1^0 x^1 x^2 ... x^n; ....;x_m^0 .... x_m^n]
	// i.e. m = 5, n = 2 x: 0, 1, 2, 3, 4
	// [0^0 0^1 0^2  
	// 1^0 1^1 1^2 
	// 2^0 2^1 2^2 
	// 3^0 3^1 3^2 
	// 4^0 4^1 4^2]
	//We only use the estimate for the middle point of the moving window for smoothing. 
	//The other rows are used only for the smoothing of the endpoint of the signal, when there are fewer data points left than the window size 2n+1.
	//Must center about m

	int m_center = (-1) * (m - 1) / 2; // want to center such that m =3 : [-1,0,1] 

	//rows
	for (int a = 0; a < m; a++) {
		//columns
		for (int b = 0; b < n + 1; b++) {
			phi(a, b) = pow(m_center, b);
		}
		m_center += 1;
	}

	cout << "phi" << phi << endl;
	//Part 5: Fill y with first m points 
		//for (int i = 0; i < m; i++) {
		//	y(i,0) = data[i];
		//}
		//cout << "y = " << y << endl;


	//Part 6 Calculate Theta
		//theta = comp_coeff(phi, y, n);
		//cout << "theta = "<< theta << endl;

	//Part 7 Loop through entire file and Filter;

		/*
		MatrixXd coeff_filter(int n, int m, MatrixXd phi, MatrixXd theta, double data) {
			//if first 2 or last 2 items: take the first two values of phi*theta
			MatrixXd y_filtered(m, 1);
			y_filtered = phi * theta;

			return y_filtered;
		}
		*/

	for (int index = 0; index <= array_length; index++) {
		cout << "current index=" << index << endl;

		//FIX THIS do centering instead of this shifting

		if (index == 0) {
			for (int data_index = 0; data_index < m; data_index++) { //(0,1,2,3,4,5)
				y(data_index, 0) = data[index + data_index];
				cout << "y[index+data_index] =" << index + data_index << endl;
			}
		}

		else if (index > (m - 1) / 2) {
			int i = 0;
			for (int data_index = -(m - 1) / 2; data_index <= (m - 1) / 2; data_index++) {
				y(i, 0) = data[index + data_index];
				i++;
			}
			cout << "else y:" << y << endl;
		}


		//cout << "y_test = " << y << endl;
		theta = comp_coeff(phi, y, n);
		//cout <<"theta = " << theta << endl;
		MatrixXd output(m, 1);
		output = coeff_filter(n, m, phi, theta);
		//first 2 elements
		if (index == 0) {
			//at (m-1)/2 element in list: 
				//write the first half values of phi*theta
			for (int x = 0; x <= (m - 1) / 2; x++) {
				cout << "if: Current index = " << index << endl << output(m - 1) / 2 << endl;
				myfile << output(x) << endl;
			}
		}
		else if ((index <= (m - 1) / 2 && index != 0) || (index > (array_length - (m - 1) / 2) && (index != (array_length - (m - 1) / 2))))
		{
			cout << "test" << endl;

		}
		//myfile << output(0)<<endl;
		//myfile << output(1) << endl;
		//myfile << output(2) << endl;
		//cout << "index = " << index << endl << "output =" << output << endl;

		else if (index == (array_length - (m - 1) / 2) - 1) {
			//at n-2 element in list 2046
				//write the last 3 values of phi*theta
			//cout << "testing else if" << endl;
			//cout << output << endl;
			for (int x = (m - 1) / 2; x < m; x++) {
				cout << "elif: Current index = " << index << endl << output(x) << endl;
				myfile << output(x) << endl;
			}

			break;
		}
		else {
			cout << "Current index = " << index << endl << output(m - 1) / 2 << endl;
			myfile << output((m - 1) / 2) << endl;
		}


		//shift y to next element

		//last 2 elements
	}
	//for (int index = 0; index < sizeof(data);index++)

	cout << "Number of elements in data= " << sizeof(data) / sizeof(double) << endl;
	cout << "array_length" << array_length << endl;

	//Close write file
	myfile.close();
	//Verification:
		//cout << data.end() << endl;
		//cout << data[2047] << endl;
		//cout << "phi*theta = " << endl << phi * theta << endl;

	cout << "Finish" << endl;


}