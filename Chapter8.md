# Introduction to Nonlinear Optimization

**Amir Beck** (https://sites.google.com/site/amirbeck314/)

- Professor of School of Mathematical Sciences, Tel-Aviv University

**References**

- **Introduction to Nonlinear Optimization, Second Edition** (2023)
  - MOS-SIAM Series on Optimization, ISBN 9781611977615 (paperback)





## Chapter 8. Convex Optimization



### 8.1 Definition

A **convex optimization problem** (or just a **convex problem**) is a problem consisting of minimizing a convex function over a convex set.
$$
\begin{align*}
\text{minimize}  &\qquad f(\mathbf{x}) \\
\text{subject to} &\qquad g_i(\mathbf{x})\leq 0, \quad 1\leq i\leq m \\
&\qquad h_j(\mathbf{x})=0, \quad 1\leq j\leq p,
\end{align*} \tag{8.2}
$$
where $f,g_i$ are <u>convex functions</u> and $h_j$​​ are <u>affine functions</u>.



Let $f\colon C\to\mathbb{R}$ be a <u>convex function</u> defined on the <u>convex set</u> $C$.

> __Theorem 8.1 (local = global in convex optimization)__ If $\mathbf{x}^*\in C$ is a <u>local minimum</u> of $f$ over $C$, then $\mathbf{x}^*$ is a <u>global minimum</u> of $f$ over $C$​.

> __Theorem 8.2__ Suppose that $f$ is <u>strictly convex</u> over $C$. If $\mathbf{x}^*\in C$ is a <u>local minimum</u> of $f$ over $C$, then $\mathbf{x}^*$ is a <u>strict (unique) global minimum</u> of $f$ over $C$.

> __Theorem 8.3 (convexity of the optimal set in convex optimization)__ The set of optimal solutions of the problem $\min\{f(\mathbf{x})\mid\mathbf{x}\in C\}$ is convex. If, in addition, $f$ is strictly convex over $C$​, then there exists at most one optimal solution.



__Example 8.4__ The problem $\min\{-2x+y\mid x^2+y^2\leq 3\}$ is convex, since $x^2+y^2\leq 3$ is a convex function. On the other hand, the problem $\min\{x^2-y\mid x^2+y^2=3\}$ is nonconvex, since $x^2+y^2=3$ is not an affine function.



### 8.2 Examples

#### 8.2.1 Linear Programming

A **linear programming (LP) problem** is an optimization problem consisting of minimizing a <u>linear objective function</u> subject to <u>linear equalities and inequalities</u>:
$$
\begin{align*}
\text{minimize} & \quad\mathbf{c}^T\mathbf{x} \\
\text{subject to} & \quad\mathbf{Ax}\leq\mathbf{b} \\
& \quad\mathbf{Bx}=\mathbf{g},
\end{align*}
$$
where $\mathbf{A}\in\mathbb{R}^{m\times n}$, $\mathbf{b}\in\mathbb{R}^{m}$, $\mathbf{B}\in\mathbb{R}^{p\times n}$, $\mathbf{g}\in\mathbb{R}^{p}$, and $\mathbf{c}\in\mathbb{R}^{n}$. This is of course a convex optimization problem.



#### 8.2.2 Convex Quadratic Problems

> __Theorem 7.10 (convexity and strict convexity of quadratic functions with positive semidefinite matrices)__ Let $f\colon\mathbb{R}^n\to\mathbb{R}$ be the quadratic function given by $f(\mathbf{x})=\mathbf{x}^T\mathbf{Ax}+2\mathbf{b}^T\mathbf{x}+c$, where $\mathbf{A}\in\mathbb{R}^{n\times n}$ is symmetric, $\mathbf{b}\in\mathbb{R}^n$, and $c\in\mathbb{R}$. Then$f$ is (resp. strictly) convex <u>if and only if</u> $\mathbf{A}\succeq\mathbf{0}$ (resp. $\mathbf{A}\succ\mathbf{0}$).

**Convex quadratic problems** are problems consisting of minimizing a <u>convex quadratic function</u> subject to <u>affine constraints</u>. A general form of problems of this class can be written as
$$
\begin{align*}
\text{minimize} & \quad\mathbf{x}^T\mathbf{Qx} + 2\mathbf{b}^T\mathbf{x} \\
\text{subject to} & \quad\mathbf{Ax}\leq\mathbf{c},
\end{align*}
$$
where $\mathbf{Q}\in\mathbb{R}^{n\times n}$ is <u>positive semidefinite</u>, $\mathbf{b}\in\mathbb{R}^n$, $\mathbf{A}\in\mathbb{R}^{m\times n}$, and $\mathbf{c}\in\mathbb{R}^m$. A well-known example of a convex quadratic problem arises in the area of linear classification.



#### 8.2.6 Convex QCQPs

A **quadratically constrained quadratic problem (QCQP)** is a problem consisting of minimizing a <u>quadratic function</u> subject to <u>quadratic inequalities and equalities</u>:
$$
\begin{align*}
\text{minimize} & \quad\mathbf{x}^T\mathbf{A}_0\mathbf{x} + 2\mathbf{b}_0^T\mathbf{x} + c_0 \\
\text{subject to} & \quad\mathbf{x}^T\mathbf{A}_i\mathbf{x} + 2\mathbf{b}_i^T\mathbf{x} + c_i \leq 0, \quad 1\leq i\leq m \\
& \quad\mathbf{x}^T\mathbf{A}_j\mathbf{x} + 2\mathbf{b}_j^T\mathbf{x} + c_j = 0, \quad m+1\leq j\leq m+p.
\end{align*}
$$
Obviously, <u>QCQPs are not necessarily convex problems</u>, but when there are <u>no equality constrains</u> ($p=0$) and all the matrices are <u>positive semidefinite</u> ($\mathbf{A}_i\succeq\mathbf{0}$ for $0\leq i\leq m$), the problem is convex and is therefore called a **convex QCQP**.
