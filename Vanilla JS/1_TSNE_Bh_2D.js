function TSNE(X, DesiredPerplexity, numberOfIterations, LearningRatio, Momentum, TradeOff){
    LearningRatio = LearningRatio * 4; //By definition of dydt ... It doesn't make sense using the *4 inside the loop
    TradeOff = TradeOff * 1.41; //sqrt of two. Relationship between de side lenght and the diagonal of the square quadtree.
    let numberOfSamplesInX = X.length;
    let numberOfDimentions = X[1].length;
    let y=zeros([numberOfSamplesInX-1,2-1]);
    let oldy=zeros([numberOfSamplesInX-1,2-1]);
    //Compute Non Symetric Pair Wise Affinities [pj|i] with perplexity [Per]
    //pj|i=exp(-||xi-xj||^2/2Sigma^2)/Sum(exp(-||xi-xk||^2/2Sigma^2))
    //Perp(Pi)=2^H(Pi)
    //H(Pi)=-Sum(pj|i)*log2(pj|i)
    //pj|i and later will be pji
    let top1 = Array(numberOfSamplesInX).fill(0);
    let top2 = Array(numberOfSamplesInX).fill(0);
    let middle1 = Array(numberOfSamplesInX).fill(0);
    let middle2 = Array(numberOfSamplesInX).fill(0);
    let bottom1 = Array(numberOfSamplesInX).fill(0);
    let bottom2 = Array(numberOfSamplesInX).fill(0);
    let p = []; //pj|i
    let IndexElements=[];
    for(let i = 0; i < X.length; i++){
        IndexElements[i] = i;
    }
    let VantagePoint = new VantagePointElement();
    VantagePoint.SelectASeedAndFindMu(X, IndexElements);
    VantagePoint.SearchKNeighbors(X, numberOfSamplesInX, DesiredPerplexity);
    let aux = getMeThePairWiseAffinities1(X, numberOfSamplesInX, numberOfDimentions);
    for(let i = 0; i <= numberOfSamplesInX - 1; i++){
        top1[i] = -50;
        bottom1[i] = -0.001;
    }
    let IShouldStay = true;
    for(let iter = 0; iter <= 50; iter++){
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
    for(let iter = 0; iter <= 50; iter++){
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
    for(let iter = 0; iter <= 100; iter++){
        //Decision Maker (you can do better than this, see it later)
        for(let i = 0; i <= numberOfSamplesInX - 1; i++){
            if(Math.abs(middle2[i] - DesiredPerplexity) < 0.01){
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
    //Set pij= ( pj|i + pi|j ) / 2 * n
    for(let i = 0; i <= numberOfSamplesInX - 2; i++){
        for(let z = 0; z <= VantagePointQueryArray[i].length - 1; z++){
            j = VantagePointQueryArray[i][z] - i - 1;
            if(j <= numberOfSamplesInX - i - 2 && j > 0){
                p[i][z] = p[i][z] / numberOfSamplesInX;
            }
        }
    }
    //Sample Initial Solution Y
    for(let i = 0; i <= numberOfSamplesInX - 2; i++){
        for(let n = 0; n <= 1; n++){
            y[i][n] = Math.random();
        }
    }
    //for t=1 to T do
    //compute low-dimensional affinities qij
    //qij = (1+||yi-yj||^2)^-1) / Sum((1+||yk-yl||^2)^-1)
    //Compute gradient dCdY
    //dCdY=4*Sum((pij-qij)*(yi-yj)*(1+||yi-yj||^2))^-1
    //Set y(t)=y(t-1) + n dCdy + a(t) * (y(t-1)-y(t-2))
    //end
    YUpload(p, y, oldy, numberOfSamplesInX, 20, Momentum * 0.5, LearningRatio, VantagePointQueryArray, TradeOff);
    YUpload(p, y, oldy, numberOfSamplesInX, numberOfIterations - 20, Momentum, LearningRatio, VantagePointQueryArray, TradeOff);
    return y;
}
function getMeThePairWiseAffinities1(X, numberOfSamplesInX, numberOfDimentions){
    //This is the part of the code that is independant of minusTwoSigmaSquared.
    //I can play this part of the code only once during the sigma search.
    let p = Array(numberOfSamplesInX).fill(0); //pj|i
    for(let i = 0; i <= numberOfSamplesInX - 2; i++){
        p[i] = Array(VantagePointQueryArray[i].length).fill(0);
        for(let z = 0; z <= VantagePointQueryArray[i].length - 1; z++){
            j = VantagePointQueryArray[i][z] - i - 1;
            if(j <= numberOfSamplesInX - i - 2 && j > 0){
                for(let n = 0; n <= numberOfDimentions - 1; n++){
                    p[i][z] = p[i][z] + Math.pow(X[i][n] - X[j + i + 1][n], 2);
                }
            }
        }
    }
    return p;
}
function getMeThePairWiseAffinities2(auxiliar, numberOfSamplesInX, minusTwoSigmaSquared){
    // get the sum of pair wise affinities and the non normilized pair wise affinities
    let sumOfPairWiseAffinities = Array(numberOfSamplesInX).fill(0);
    let p = Array(numberOfSamplesInX).fill(0); //pj|i
    for(let i = 0; i <= numberOfSamplesInX - 2; i++){
        p[i] = Array(VantagePointQueryArray[i].length).fill(0);
        for(let z = 0; z <= VantagePointQueryArray[i].length - 1; z++){
            j = VantagePointQueryArray[i][z] - i - 1;
            if(j <= numberOfSamplesInX - i - 2 && j > 0){
                p[i][z] = Math.exp(auxiliar[i][z] / (minusTwoSigmaSquared[i] + 0.000001));
                sumOfPairWiseAffinities[i] = sumOfPairWiseAffinities[i] + p[i][z];
                sumOfPairWiseAffinities[j + i + 1] = sumOfPairWiseAffinities[j + i + 1] + p[i][z];
            }
        }
    }
    // get the normilized pair wise affinities
    for(let i = 0; i <= numberOfSamplesInX - 2; i++){
        for(let z = 0; z <= VantagePointQueryArray[i].length - 1; z++){
            j = VantagePointQueryArray[i][z] - i - 1;
            if(j <= numberOfSamplesInX - i - 2 && j > 0){
                p[i][z] = p[i][z] / (sumOfPairWiseAffinities[j + i + 1] + 0.000001);
            }
        }
    }
    return p;
}
function GetMeThePerplexity(p, numberOfSamplesInX){
    let Perplexities=Array(numberOfSamplesInX).fill(0);
    for(i = 0; i <= numberOfSamplesInX - 2; i++){
        for(let z = 0; z <= VantagePointQueryArray[i].length - 1; z++){
            j = VantagePointQueryArray[i][z] - i - 1;
            if(j <= numberOfSamplesInX - i - 2 && j > 0){
                let aux = p[i][z] * Math.log2(p[i][z] + 0.001);
                Perplexities[i] = Perplexities[i] + aux;
                Perplexities[j + i + 1] = Perplexities[j + i + 1] + aux;
            }
        }
    }
    for(i = 0; i <= numberOfSamplesInX - 1; i++){
        Perplexities[i] = Math.pow( 2, -Perplexities[i]);
    }
    return Perplexities;
}
function YUpload(p, y, oldy, numberOfSamplesInX, numberOfIterations, Momentum, LearningRatio, VantagePointQueryArray, TradeOff){
    let aux;
    let aux1;
    let j;
    let BiggestY = 2;
    let IndexElements=[];
    let ResultQT = new QuadTreeResults();
    for(let i = 0; i < y.length; i++){
        IndexElements[i] = i;
    }
    for(let iter = 0; iter <= numberOfIterations; iter++){
        let Fattr = zeros([numberOfSamplesInX-1,1]);
        let Frep = zeros([numberOfSamplesInX-1,1]);
        let Sumq = 0;
        QuadTree = new QuadtreeElement([0,0], BiggestY);
        QuadTree.InsertInBoxes(y, IndexElements);
        BiggestY = 0;
        //F attr
        for(let i = 0; i <= numberOfSamplesInX - 2; i++){
            for(let z = 0; z <= VantagePointQueryArray[i].length - 1; z++){
                j = VantagePointQueryArray[i][z] - i - 1
                if(j <= numberOfSamplesInX - i - 2 && j > 0){
                    aux1=p[i][z] * CalculateZQij( i, j + i + 1);
                    for(let n = 0; n <= 1; n++){
                        aux = aux1 * (y[i][n] - y[j+ i + 1][n]);
                        Fattr[i][n] = Fattr[i][n] + aux;
                        Fattr[j + i + 1][n] = Fattr[j + i + 1][n] - aux;
                    }
                }
            }
        }
        for(let i = 0; i <= numberOfSamplesInX - 2; i++){
            QuadTree.ListOfEquivalentBodiesOfI(y, i, TradeOff, ResultQT);
            for(let z = 0; z <= ResultQT.ResultOfTheQueryQT1.length - 1; z++){
                aux = CalculateZQij2( i, ResultQT.ResultOfTheQueryQT1[z]);
                Sumq = Sumq + ResultOfTheQueryQT2[z] * 2 / aux;
                aux1 = Math.pow(aux, 2);
                for(let n = 0; n <= 1; n++){
                    aux = aux1 * (y[i][n] - ResultQT.ResultOfTheQueryQT1[z][n] * ResultOfTheQueryQT2[z])
                    Frep[i][n] = Frep[i][n] - aux;
                }
            }
            for(let z = 0; z <= ResultQT.ResultOfTheQueryQT3.length - 1; z++){
                aux = CalculateZQij( i, ResultQT.ResultOfTheQueryQT3[z]);
                Sumq = Sumq +  2 / aux;
                aux1 = Math.pow(aux, 2);
                for(let n = 0; n <= 1; n++){
                    aux = aux1 * (y[i][n] - y[ResultQT.ResultOfTheQueryQT3[z]][n]);
                    Frep[i][n] = Frep[i][n] - aux;
                }
            }
        }
        //y adjustment
        for(let i = 0; i <= numberOfSamplesInX - 1; i++){
            for(let n = 0; n <= 1; n++){
                aux = y[i][n];
                y[i][n] = y[i][n] - LearningRatio * (Fattr[i][n] + Frep[i][n] / Sumq) + Momentum * (y[i][n] - oldy[i][n]);
                oldy[i][n] = aux;
                if(BiggestY < Math.abs(y[i][n])){BiggestY = Math.abs(y[i][n])}
            }
        }
    }
    function CalculateZQij( i, j){
        let aux=0;
        for(let n = 0; n <= 1; n++){
            aux = aux + Math.pow(y[i][n] - y[j][n], 2);
        }
        return Math.pow(1 + aux, -1);
    }
    function CalculateZQij2( i, CenterOfMass){
        let aux=0;
        for(let n = 0; n <= 1; n++){
            aux = aux + Math.pow(y[i][n] - CenterOfMass[n], 2);
        }
        return Math.pow(1 + aux, -1);
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
