let DesiredPerplexity = 6;
let numberOfIterations = 1;
let LearningRatio = 1;
let Momentum = 0.8;
let p = [];
let y = [];
let oldy = [];
let numberOfSamplesInX;
let Maximum = -100;
let Minimum = 100;
let RangeX = 0;
let numberOfDimentionsInLowDimensionalSpace = 2;
function setup() {
    createCanvas(400, 400);
    LearningRatio = LearningRatio * 4; //By definition of dydt ... It doesn't make sense using the *4 inside the loop
    numberOfSamplesInX = X.length;
    let numberOfDimentions = X[1].length;
    y=zeros([numberOfSamplesInX-1,numberOfDimentionsInLowDimensionalSpace-1]);
    oldy=zeros([numberOfSamplesInX-1,numberOfDimentionsInLowDimensionalSpace-1]);
    //Compute Non Symetric Pair Wise Affinities [pj|i] with perplexity [Per]
    //pj|i=exp(-||xi-xj||^2/2Sigma^2)/Sum(exp(-||xi-xk||^2/2Sigma^2))
    //Perp(Pi)=2^H(Pi)
    //H(Pi)=-Sum(pj|i)*log2(pj|i)
    //pj|i and later will be pji
    p = SearchMeForFixedPerpexity(X, numberOfSamplesInX, numberOfDimentions, DesiredPerplexity);
    //Set pij= ( pj|i + pi|j ) / 2 * n
    for(let i = 0; i <= numberOfSamplesInX - 2; i++){
        for(let j = 0; j <= numberOfSamplesInX - i - 2; j++){
            p[i][j] = p[i][j] / numberOfSamplesInX;
        }
    }
    //Sample Initial Solution Y
    for(let i = 0; i <= numberOfSamplesInX - 2; i++){
        for(let n = 0; n <= numberOfDimentionsInLowDimensionalSpace - 1; n++){
            y[i][n] = Math.random()-0.5;
        }
    }
    //for t=1 to T do
    //compute low-dimensional affinities qij
    //qij = (1+||yi-yj||^2)^-1) / Sum((1+||yk-yl||^2)^-1)
    //Compute gradient dCdY
    //dCdY=4*Sum((pij-qij)*(yi-yj)*(1+||yi-yj||^2))^-1
    //Set y(t)=y(t-1) + n dCdy + a(t) * (y(t-1)-y(t-2))
    //end
}
function draw() {
    background(0);
    YUpload(p, y, oldy, numberOfSamplesInX, numberOfDimentionsInLowDimensionalSpace, numberOfIterations, Momentum, LearningRatio);
    Maximum = -100;
    Minimum = 100;
    for(let i = 0; i < X.length; i++){
        if(y[i][0]<Minimum){Minimum = y[i][0]}
        if(y[i][1]<Minimum){Minimum = y[i][1]}
        if(y[i][0]>Maximum){Maximum = y[i][0]}
        if(y[i][1]>Maximum){Maximum = y[i][1]}
    }
    RangeX = Maximum - Minimum;
    for(let i = 0; i < y.length; i++){
        stroke('purple'); // Change the color
        strokeWeight(3); // Make the points 10 pixels in      
        point(y[i][0]*300/RangeX+200, y[i][1]*300/RangeX+200);   
    }
}
function getMeThePairWiseAffinities1(X, numberOfSamplesInX, numberOfDimentions){
    //This is the part of the code that is independant of minusTwoSigmaSquared.
    //I can play this part of the code only once during the sigma search.
    let p = GenerateHalfAMatrix(numberOfSamplesInX); //pj|i
    for(let i = 0; i <= numberOfSamplesInX - 2; i++){
        for(let j = 0; j <= numberOfSamplesInX - i - 2; j++){
            for(let n = 0; n <= numberOfDimentions - 1; n++){
                p[i][j] = p[i][j] + Math.pow(X[i][n] - X[j + i + 1][n], 2);
            }
        }
    }
    return p;
}
function getMeThePairWiseAffinities2(auxiliar, numberOfSamplesInX, minusTwoSigmaSquared){
    // get the sum of pair wise affinities and the non normilized pair wise affinities
    let sumOfPairWiseAffinities = Array(numberOfSamplesInX).fill(0);
    let p = GenerateHalfAMatrix(numberOfSamplesInX);
    for(let i = 0; i <= numberOfSamplesInX - 2; i++){
        for(let j = 0; j <= numberOfSamplesInX - i - 2; j++){
            p[i][j] = Math.exp(auxiliar[i][j] / (minusTwoSigmaSquared[i] + 0.000001));
            sumOfPairWiseAffinities[i] = sumOfPairWiseAffinities[i] + p[i][j];
            sumOfPairWiseAffinities[j + i + 1] = sumOfPairWiseAffinities[j + i + 1] + p[i][j];
        }
    }
    // get the normilized pair wise affinities
    for(let i = 0; i <= numberOfSamplesInX - 2; i++){
        for(let j = 0; j <= numberOfSamplesInX - i - 2; j++){
            p[i][j] = p[i][j] / (sumOfPairWiseAffinities[j + i + 1] + 0.000001);
        }
    }
    return p;
}
function SearchMeForFixedPerpexity(X, numberOfSamplesInX, numberOfDimentions, DesiredPerplexity){
    let top1 = Array(numberOfSamplesInX).fill(0);
    let top2 = Array(numberOfSamplesInX).fill(0);
    let middle1 = Array(numberOfSamplesInX).fill(0);
    let middle2 = Array(numberOfSamplesInX).fill(0);
    let bottom1 = Array(numberOfSamplesInX).fill(0);
    let bottom2 = Array(numberOfSamplesInX).fill(0);
    let p = []; //pj|i
    let aux = getMeThePairWiseAffinities1(X, numberOfSamplesInX, numberOfDimentions);
    for(let i = 0; i <= numberOfSamplesInX - 1; i++){
        top1[i] = -50;
        bottom1[i] = -0.001;
    }
    let IShouldStay = true;
    for(let iter = 1; iter <= 50; iter++){
        IShouldStay = true;
        p = getMeThePairWiseAffinities2(aux, numberOfSamplesInX, top1);
        top2 = GetMeThePerplexity(p, numberOfSamplesInX);
        for(let i = 0; i <= numberOfSamplesInX - 1; i++){
            if(top2[i] < DesiredPerplexity){
                IShouldStay = false;
                top1[i] = top1[i] * 2;
            }
        }
        if(IShouldStay){iter = 51};
    }
    for(let iter = 1; iter <= 50; iter++){
        IShouldStay = true;
        p = getMeThePairWiseAffinities2(aux, numberOfSamplesInX, bottom1);
        bottom2 = GetMeThePerplexity(p, numberOfSamplesInX);
        for(let i = 0; i <= numberOfSamplesInX - 1; i++){
            if(bottom2[i] > DesiredPerplexity){
                IShouldStay = false;
                bottom1[i] = bottom1[i] * 0.5;
            }
        }
        if(IShouldStay){iter = 51};
    }
    for(let i = 0; i <= numberOfSamplesInX - 1; i++){
        middle1[i] = (top1[i] + bottom1[i]) / 2;
    }
    p = getMeThePairWiseAffinities2(aux, numberOfSamplesInX, middle1);
    middle2 = GetMeThePerplexity(p, numberOfSamplesInX);
    for(let iter = 1; iter <= 100; iter++){
        //Decision Maker (you can do better than this, see it later)
        for(let i = 0; i <= numberOfSamplesInX - 1; i++){
            if(Math.abs(middle2[i] - DesiredPerplexity) < 0.01 || iter == 500){
            }else if(middle2[i] > DesiredPerplexity){
                top1[i] = middle1[i];
                top2[i] = middle2[i];
                middle1[i] = (top1[i] + bottom1[i]) / 2;
            }else{
                bottom1[i] = middle1[i];
                bottom2[i] = middle2[i];
                middle1[i] = (top1[i] + bottom1[i]) / 2;
            }
        }
        p = getMeThePairWiseAffinities2(aux, numberOfSamplesInX, middle1);
        middle2 = GetMeThePerplexity(p, numberOfSamplesInX);
    }
    return p;
}
function GetMeThePerplexity(p, numberOfSamplesInX){
    let Perplexities=Array(numberOfSamplesInX).fill(0);
    for(i = 0; i <= numberOfSamplesInX - 2; i++){
        for(j = 0; j <= numberOfSamplesInX - i - 2; j++){
            let aux = p[i][j] * Math.log2(p[i][j] + 0.001);
            Perplexities[i] = Perplexities[i] + aux;
            Perplexities[j + i + 1] = Perplexities[j + i + 1] + aux;
        }
    }
    for(i = 0; i <= numberOfSamplesInX - 1; i++){
        Perplexities[i] = Math.pow( 2, -Perplexities[i]);
    }
    return Perplexities;
}
function YUpload(p, y, oldy, numberOfSamplesInX, numberOfDimentionsInLowDimensionalSpace, numberOfIterations, Momentum, LearningRatio){
    let aux;
    let aux1;
    for(let iter = 0; iter <= numberOfIterations; iter++){
        let Zq = GenerateHalfAMatrix(numberOfSamplesInX);
        let Sumq = 0;
        let dCdYi = zeros([numberOfSamplesInX-1,numberOfDimentionsInLowDimensionalSpace-1]);
        for(let i = 0; i <= numberOfSamplesInX - 2; i++){
            for(let j = 0; j <= numberOfSamplesInX - i - 2; j++){
                for(let n = 0; n <= numberOfDimentionsInLowDimensionalSpace - 1; n++){
                    Zq[i][j] = Zq[i][j] + Math.pow(y[i][n] - y[j + i + 1][n], 2);
                }
                Zq[i][j] = Math.pow(1 + Zq[i][j], -1);
                Sumq = Sumq + 2 * Zq[i][j];
            }
        }
        for(let i = 0; i <= numberOfSamplesInX - 2; i++){
            for(let j = 0; j <= numberOfSamplesInX - i - 2; j++){
                aux1 = (p[i][j] - Zq[i][j] / Sumq) * Zq[i][j];
                for(let n = 0; n <= numberOfDimentionsInLowDimensionalSpace - 1; n++){
                    aux = aux1 * (y[i][n] - y[j+ i + 1][n]);
                    dCdYi[i][n] = dCdYi[i][n] + aux;
                    dCdYi[j + i + 1][n] = dCdYi[j + i + 1][n] - aux;
                }
            }
        }
        //y adjustment
        for(let i = 0; i <= numberOfSamplesInX - 1; i++){
            for(let n = 0; n <= numberOfDimentionsInLowDimensionalSpace - 1; n++){
                aux = y[i][n];
                y[i][n] = y[i][n] - LearningRatio * dCdYi[i][n] + Momentum * (y[i][n] - oldy[i][n]);
                oldy[i][n] = aux;
            }
        }
    }
}
function GenerateHalfAMatrix(dimensions){
    let aux = Array(dimensions - 1).fill(0);
    let aux1 = [];
    for(let i=0; i <= dimensions - 2; i++){
        aux1=Array(dimensions - i - 1).fill(0);
        aux[i]=aux1;
    }
    return aux;
}
function zeros(dimensions){
    var array = [];
    for (var i = 0; i <= dimensions[0]; ++i) {
        array.push(dimensions.length == 1 ? 0 : zeros(dimensions.slice(1)));
    }
    return array;
}