#ifndef ELEMENTO_H
	#define ELEMENTO_H
	/*
	 * elemento.cpp
	 *
	 *  Created on: Mar 23, 2011
	 *      Author: oscar
	 */

//	using namespace std;

	namespace Tera{
	/**
	*	Función para calcular las fuerzas en un elemento triangular
	*   @param double referencia[3][3], Vertices iniciales del elemento
	*	@param double deformado[3][3], Vertices del elemento deformado en el plano
	*	@return double** fuerzas, Apuntador hacia las fuerzas calculadas en cada nodo
	*/

	void fuerzas(double referencia[3][3], double deformado[3][3], double fuerzas[3][3], double ks);


	/**
	*	Función para realizar cross product entre dos vectores
	*   @param double a[3], Vector a
	*	@param double b[3], Vector b
	*	@param double &resultado[3], Paso por parametros
	*/
	void cross(double a[3], double b[3], double &x, double &y, double &z);

	/**
	*	Función para calcular la norma de un vector
	*   @param double a[3], Vector a
	*	@return double norma, Norma del vector
	*/
	double norm(double a[3]);

	/**
	*	Función para calcular el producto Matriz dot Vector
	*   @param M[3][3], Matriz de tamaño 3x3
	*   @param V[3], Vector de tamaño 3
	*   @param V[3], Vector de tamaño 3
	*/
	void MdotV(double M[3][3], double V[3], double &x, double &y, double &z);

	/**
	*	Función para llevar los elementos de referencia y deformado a una misma base
	*   @param double referencia[3][3], Vertices iniciales del elemento
	*	@param double deformado[3][3], Vertices del elemento deformado en 3D
	*	@return double** vertices, Apuntador hacia los vertices en un plano común
	*/

	void rotacion(double referencia[3][3], double deformado[3][3], double nfuerzas[3][3], double ks);
	

}//end namespace
#endif
