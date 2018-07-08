// stacks stuff.cpp : Defines the entry point for the console application.
//
//


//TODO: lots of error catching required
//TODO: add more comments
//TODO: tweak to allow exponentials
//TODO: auto add '*'s
//TODO: the rest of it
//TODO: add function getLimits, which asks for input upper and lower bounds as well as checking for e's and pi's


#include "stdafx.h"
#include <iostream>
#include <string>
#include <string.h>
#include <stack>
#include <queue>
#include <sstream>
#include <cctype>
#include <cmath>
#include <cstdlib>

using namespace std;

bool isFunction(string function) //checks if string is trig function
{
	cout << "isFunction started\n";
	if (function == "sin(" || function == "cos(" || function == "tan(") return true;
	else return false;
}

string queueToString(queue<string> myQueue) //converts queue to string, delimited by semi-colons
{
	cout << "queueToString started\n";
	string myString = "";

	while (!myQueue.empty()) //appends each item of the queue to string with semi-colon in between
	{
		myString += myQueue.front();
		myString += "; ";
		myQueue.pop();
	}
	
	return myString;
}

string intToString(signed long long int num) //converts an integer to a string to be stored by the program
{
	cout << "intToString started\n";
	ostringstream convert; 
	convert << num; //stringstream set up to which the integer is sent to
	return convert.str(); //stringstream converted to string and returned
}

int getRanking(string operation) /*gets the relative rankings of different operations, 
higher number meaning higher precedence*/
{
	cout << "getRanking started\n";
	if (operation == "+" || operation == "-")
	{
		return 1;
	}
	else if (operation == "/" || operation == "*")
	{
		return 2;
	}
	else if (operation == "(") //open bracket always has lowest precedence
	{
		return 0;
	}
	else if (operation == "^")
	{
		return 4;
	}
	else
	{
		return -1; //error return code, operation should be one of the above
	}
}

queue<string> makePostfix(queue<string> inQueue) //converts a queue in infix notation to a queue in postfix notation
{
	cout << "makePostfix started\n";
	queue<string> outQueue;
	stack<string> operators; //stack to contain all the operators in the expression


	cout << queueToString(inQueue) << endl;
	
	while (!inQueue.empty()) //looping through the queue to check each value
	{
		if (inQueue.front() == "(")
		{
			operators.push(inQueue.front());
		}
		else if (inQueue.front() == ")")
		{
			while (operators.top() != "(" && isFunction(operators.top()) == false) /*if a close bracket is found, all 
			operators are popped to the output until the corresponding open bracket is found*/
			{
				outQueue.push(operators.top());
				operators.pop();
			}
			if (isFunction(operators.top()))
			{
				outQueue.push(operators.top()); 
				//outQueue.push("**");
			}
			operators.pop();
		}
		else if (inQueue.front() == "+" || inQueue.front() == "-" || inQueue.front() == "/" || inQueue.front() == "*")
		{
			if (operators.empty() == true || operators.top() == "(") /*if there are no operators before this, 
				or if it is part of a pair of brackets*/
			{
				operators.push(inQueue.front());
			}
			else
			{
				while (getRanking(operators.top()) >= getRanking(inQueue.front())) /* if there is an operator in the stack 
					with higher precedence than this one*/
				{
					outQueue.push(operators.top());
					operators.pop(); //top of stack pushed to output
					if (operators.empty() == true) /*if the stack is empty, it means that the expression has not been input
						properly, as there wouldn't be the right number of numbers*/
					{
						break; //TODO: fix error catching here
					}
				}
				operators.push(inQueue.front()); /*once all higher precedence operators are gone from the stack, the operator 
				is pushed to it*/
			}
		}
		else if (inQueue.front() == "^")
		{
			if (operators.empty() == true || operators.top() == "(") /*if there are no operators before this,
																	 or if it is part of a pair of brackets*/
			{
				operators.push(inQueue.front());
			}
			else
			{
				while (getRanking(operators.top()) > getRanking(inQueue.front())) /* if there is an operator in the stack 
																with higher precedence than (but not equal to) this one*/
				{
					outQueue.push(operators.top());
					operators.pop(); //top of stack pushed to output
					if (operators.empty() == true) /*if the stack is empty, it means that the expression has not been input
												   properly, as there wouldn't be the right number of numbers*/
					{
						break; //TODO: fix error catching here
					}
				}
				operators.push(inQueue.front()); /*once all higher precedence operators are gone from the stack, the operator
												 is pushed to it*/
			}
		}
		else if (inQueue.front() == "sin(" || inQueue.front() == "tan(" || inQueue.front() == "cos(")
		{
			operators.push(inQueue.front());
		}
		else //TODO: add error checking here
		{
			outQueue.push(inQueue.front()); //if any other (valid) character, push to output
		}

		inQueue.pop();
	} 

	while (operators.empty() == false) //while there are still items on the stack, push them to the output
	{
		outQueue.push(operators.top());
		operators.pop();
	}

	cout << queueToString(outQueue) << endl; //debug print of resulting queue 

	return outQueue;
}

queue<string> makeQueue() //converts the input string into a queue
{
	cout << "makeQueue started\n";

	queue<string> tokens; 
	string instring;
	int polarity = 1;

	cout << "Enter Expression:\n";
	getline(cin, instring);

	signed long long int currentNumber = 0; //this needs to be large to store very long numbers
	string currentToken = "";
	string currentString = ""; //for trig  functions to be stored in

	for (int i = 0; i < instring.length(); ++i)
	{
		currentToken = instring.substr(i, 1);
		if (currentToken == "+" || currentToken == "-" || currentToken == "/" ||currentToken == "*" ||
			currentToken == "(" || currentToken == ")" || currentToken == "^" || currentToken == "x"|| currentToken == "e")
		{
			if (currentToken != "-" && currentNumber == 0)
			{
				if (currentToken == "(" && currentString != "")
				{
					tokens.push(currentString + currentToken);
					currentString = "";
				}
				else
				{
					tokens.push(currentToken);
				}
			}
			if (currentNumber == 0 && currentToken == "-" && instring.substr(i-1,1) != ")") polarity = -1;
			else if (currentNumber != 0 ) //when a operator is reached, the current number is pushed to the queue
			{
				tokens.push(intToString(polarity * currentNumber));
				polarity = 1;
				currentNumber = 0;
				tokens.push(currentToken);
			}
			else if (currentToken == "-" && instring.substr(i - 1, 1) == ")")
			{
				tokens.push(currentToken);
			}
			
		}
		else if (isdigit(currentToken.c_str()[0]))
		{
			currentNumber = 10 * currentNumber + atoi(currentToken.c_str()); /*if character is another digit, 
			add it to the current number*/
		}
		else if (currentToken == "s" || currentToken == "i" || currentToken == "n" || currentToken == "c" || 
			currentToken == "o" || currentToken == "t" || currentToken == "a"|| currentToken == "p")
		{
			currentString += currentToken;
		}

		if (currentString == "pi")
		{
			tokens.push("pi");
			currentString = "";
		}

	}

	if (currentNumber != 0)
	{
		tokens.push(intToString(polarity * currentNumber));
		currentNumber = 0;
	}

	return tokens;
}

double binaryEval(double num1, double num2, string operation) /*this subroutine takes two operands and an operator
and returns the result of the calculation*/
{
	cout << "binaryEval started\n";
	if (operation == "+")
	{
		return num1 + num2;
	}
	else if (operation == "-")
	{
		return num2 - num1;
	}
	else if (operation == "*")
	{
		return num1*num2;
	}
	else if (operation == "/")
	{
		return num2 / num1;
	}
	else if (operation == "^")
	{
		return pow(num2, num1);
	}
}

double functionEval(double num, string operation)
{
	cout << "functionEval started\n";
	if (operation == "sin(")
	{
		return sin(num);
	}
	else if (operation == "cos(")
	{
		return cos(num);
	}
	else if (operation == "tan(")
	{
		return tan(num);
	}
}

double evaluateExpression(double x, queue<string> expression) /*this subroutine evaluates the value of an expression
															using the postfix queue and the value x entered*/
{
	cout << "evaluateExpression started\n";
	stack<double> operands;
	string currentToken = "";
	double num1 = 0;
	double num2 = 0;

	while (!expression.empty())
	{
		currentToken = expression.front();
		expression.pop();

		if (isdigit(currentToken.c_str()[0]) || isdigit(currentToken.c_str()[1])) //if number, it is pushed to the operands stack
		{
			operands.push(atoi(currentToken.c_str())); 
		}
		else if (currentToken == "+" || currentToken == "-" || currentToken == "/"
			|| currentToken == "*" || currentToken == "^")
		{
			num1 = operands.top();
			operands.pop();
			num2 = operands.top();
			operands.pop(); //the two numbers at the top are popped when an operator is reached

			operands.push(binaryEval(num1, num2, currentToken)); //value of two numbers combined with the operator
		}
		else if (currentToken == "x") //if the string is x, its value is pushed to the stack
		{
			operands.push(x);
		}
		else if (currentToken == "e")
		{
			operands.push(2.71828182845904523536);
		}
		else if (currentToken == "pi")
		{
			operands.push(3.14159265358979323846);
		}
		else if (currentToken == "sin(" || currentToken == "cos(" || currentToken == "tan(")
		{
			num1 = operands.top();
			operands.pop();
			operands.push(functionEval(num1, currentToken));
		}

	}

	if (operands.size() == 1) //there should only be one item on the result stack at the end, otherwise there has been an error
	{
		cout << operands.top() << endl;
		return operands.top();
	}
	else
	{
		return -1;
	}
}

int valueExpression()
{
	cout << "valueExpression started\n";
	queue<string> expression;
	int x;

	expression = makePostfix(makeQueue());

	cout << "x = ";

	cin >> x;

	evaluateExpression(x, expression);

	return 0;
}

double trapeziumRule()
{
	cout << "trapeziumRule started\n";
	queue<string> expression;
	double upper;
	double lower;
	int noIntervals;

	cout << "upper = ";
	cin >> upper;
	cout << "lower = ";
	cin >> lower;
	cout << "noIntervals = ";
	cin >> noIntervals;
	cin.ignore(256, '\n');
	expression = makePostfix(makeQueue());
	double total = 0;
	double x = lower;

	total += (evaluateExpression(lower, expression) + evaluateExpression(upper, expression));

	double h = (upper - lower) / noIntervals;

	while (!(x + h >= upper))
	{
		x = x + h;
		total = total + 2*evaluateExpression(x, expression);
	} 

	total = total * (h / 2);

	cout << total;
	return total;
}

double simpsonsRule()
{
	cout << "simpsonsRule started\n";
	queue<string> expression;
	double upper;
	double lower;
	int noIntervals;

	cout << "upper = ";
	cin >> upper;
	cout << "lower = ";
	cin >> lower;
	cout << "noIntervals = ";
	cin >> noIntervals;
	cin.ignore(256, '\n');
	expression = makePostfix(makeQueue());

	double total = 0;
	double x = lower;
	double h = (upper - lower) / noIntervals;
	bool parity = true;

	total += (evaluateExpression(lower, expression) + evaluateExpression(upper, expression));
	while (x + h < upper)
	{
		x += h;
		if (parity == true)	total += (4 * evaluateExpression(x, expression));
		else total += 2 * evaluateExpression(x, expression);

		parity = !parity;
	}

	/*if (parity == true) 
	{ 
		return -1; 
	}
	else*/
	{ 
		total *= (h / 3);
		cout << total;
		return total;
	}

}

int main()
{
	cout << "main started\n";
	//makePostfix(makeQueue());
	valueExpression();
	//trapeziumRule();
	//simpsonsRule();
	cin.ignore(256, '\n');
	cin.get();
	return 0;
}

