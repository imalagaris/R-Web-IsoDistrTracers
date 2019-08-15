___
**Variance**  
$$Var(X) = E[(X - \mu_x)^2]$$
$$Var(X) = \sum_{x_k\in\mathbb{R}}(x_k - \mu_X)^2P_X(x_k)$$
$$Var(X) = E[X^2] - [EX]^2$$
$$EX^2 = \sum_{x_k\in\mathbb{R}}x_k^2P_X(x_k)$$

*Variance is not linear*  
$$Var(\alpha X + \beta) = \alpha^2Var(X)$$  
If $X_1, X_2,..., X_n$ are **independent** random variables, then the variance of the random variable $X =  X_1, X_2,..., X_n$ is
$$Var(X) = Var(X_1) + Var(X_2) +...+ Var(X_n)$$
<br> </br>

___
**Marginal PDFs**  
$$f_X(x) = \int_{-\infty}^{\infty}f_{XY}(x,y)dy \text{ for all x,}$$
$$f_Y(y) = \int_{-\infty}^{\infty}f_{XY}(x,y)dx \text{ for all y.}$$
$$F_{XY}(x,y) = \int_{-\infty}^y\int_{-\infty}^xf_{XY}(u,v)dudv$$
$$f_{XY}(x,y) = \frac{\partial^2} {\partial x\partial y}F_{XY}(x,y)$$
<br> </br>

___
**Law of Total Probability**  
$$P(A) = \int_{-\infty}^{\infty}P(A|X = x) f_X(x) dx$$
<br> </br>

___
**Law of Total Expectation**  
$$EY = \int_{-\infty}^{\infty}E[Y|X = x]f_X(x)dx$$
$$EY = E [E [Y|X] ]$$
**Law of Total Variance**
$$Var(Y) = E[Var(Y|X)] + Var(E[Y|X])$$
$$Var(Y) \geq E[Var(X|Y)]$$
<br> </br>

___
**Convolution integral**  
If X and Y are *independent random* variables and Z = X + Y then
$$f_Z(z) = \int_{-\infty}^{\infty}f_X(w)f_Y(z-w)dw$$
<br> </br>

___
**Theorem**  
If X and Y are independent and  
$X \sim N(\mu_X,\sigma_X^2)$ and  
$Y \sim N(\mu_Y,\sigma_Y^2)$ then  
$$X + Y \sim N(\mu_X + \mu_Y, \sigma_X^2 + \sigma_Y^2)$$
<br> </br>
