let VantagePointElementCapacity = 1;
let VantagePointMaximumInception =16;
let VantagePointQueryArray = Array(31).fill(0);
let ResultOfTheQueryVP = [];
let numberOfDimentions;
class VantagePointElement{
    constructor(){
        this.IndexOfSeed = undefined;
        this.Mu = undefined;
        this.ListOfIndexElementsInsideMe = [];
        this.ListOfIndexElementsOutsideMe = [];
        this.InsideVantagePointSon = undefined;
        this.OutsideVantagePointSon = undefined;
        this.AmountOfElementsInMe = undefined;
        this.DistanceToTheSample = undefined;
        this.level = 0;
    }
    SearchKNeighbors(DataBase, numberOfSamplesInX, DesiredPerplexity){
        let threshold;
        let RequestedNeighbors = DesiredPerplexity * 3;
        for(let i = 0 ; i < numberOfSamplesInX; i++){
            threshold = 0.0001;
            for(let iter=0; iter < 100; iter++){
                ResultOfTheQueryVP = [];
                this.SearchNeighborOfI(DataBase, i, threshold);
                if (ResultOfTheQueryVP.length < RequestedNeighbors){
                    threshold = (threshold + 0.001) * 2
                }else{
                    iter = 100;
                }
            }
            VantagePointQueryArray[i] = ResultOfTheQueryVP;
        }
    }
    SearchAllNeighbors(DataBase, numberOfSamplesInX, threshold){
        for(let i = 0 ; i < numberOfSamplesInX; i++){
            ResultOfTheQueryVP = [];
            this.SearchNeighborOfI(DataBase, i, -threshold[i]);
            VantagePointQueryArray[i] = ResultOfTheQueryVP;
        }
    }
    SearchNeighborOfI(DataBase, Sample, threshold){
        //Measure my distance
        this.DistanceToTheSample = 0;
        for(let i=0; i < numberOfDimentions; i++){
            this.DistanceToTheSample = this.DistanceToTheSample + Math.pow(DataBase[this.IndexOfSeed][i]-DataBase[Sample][i],2);
        }
        this.DistanceToTheSample = Math.pow(this.DistanceToTheSample,0.5);
        //If i am far away
        if(this.DistanceToTheSample > this.Mu + threshold){
            if(typeof this.OutsideVantagePointSon != 'undefined'){
                this.OutsideVantagePointSon.SearchNeighborOfI(DataBase, Sample, threshold);
            }else{
                for(let i=0; i < this.ListOfIndexElementsOutsideMe.length; i++){
                    ResultOfTheQueryVP.push(this.ListOfIndexElementsOutsideMe[i]);
                }
            }
        }else if(this.Mu > threshold + this.DistanceToTheSample){
        //If i am very very close
            if(typeof this.InsideVantagePointSon != 'undefined'){
                this.InsideVantagePointSon.SearchNeighborOfI(DataBase, Sample, threshold);
            }else{
                for(let i=0; i < this.ListOfIndexElementsInsideMe.length; i++){
                    ResultOfTheQueryVP.push(this.ListOfIndexElementsInsideMe[i]);
                }
            }
        }else{
        //Is in between
            if(typeof this.OutsideVantagePointSon != 'undefined'){
                this.OutsideVantagePointSon.SearchNeighborOfI(DataBase, Sample, threshold);
            }else{
                for(let i=0; i < this.ListOfIndexElementsOutsideMe.length; i++){
                    ResultOfTheQueryVP.push(this.ListOfIndexElementsOutsideMe[i]);
                }
            }
            if(typeof this.InsideVantagePointSon != 'undefined'){
                this.InsideVantagePointSon.SearchNeighborOfI(DataBase, Sample, threshold);
            }else{
                for(let i=0; i < this.ListOfIndexElementsInsideMe.length; i++){
                    ResultOfTheQueryVP.push(this.ListOfIndexElementsInsideMe[i]);
                }
            }
        }
        this.DistanceToTheSample = undefined;
    }
    SelectASeedAndFindMu(DataBase, IndexOfElementsToEvaluate){
        numberOfDimentions = DataBase[0].length;
        //Select a Seed
        let TheChoosenIndex = Math.floor(Math.random()*(IndexOfElementsToEvaluate.length));
        this.IndexOfSeed = IndexOfElementsToEvaluate[TheChoosenIndex];
        this.AmountOfElementsInMe = IndexOfElementsToEvaluate.length;
        //FindTheMu
        let NonDefinitiveListOfIndexInsideMe=[];
        let NonDefinitiveListOfIndexOutsideMe=[];
        let TargetOfElementsInsideMe = Math.floor(this.AmountOfElementsInMe/2);
        let aux;
        let TopMu = 400;
        let BottomMu = 0;
        let CounterOfElementsInsideMe = 0;
        let ToleranceOfIteration = 100;
        let MiddleMu = 10;
        for(let iter = 0; iter < ToleranceOfIteration; iter++){
            CounterOfElementsInsideMe = 0;
            aux = Array(this.AmountOfElementsInMe).fill(0);
            for(let Sample = 0; Sample < this.AmountOfElementsInMe; Sample++){
                for(let i=0; i < numberOfDimentions; i++){
                    aux[Sample] = aux[Sample] + Math.pow(DataBase[this.IndexOfSeed][i]-DataBase[IndexOfElementsToEvaluate[Sample]][i],2);
                }
                aux[Sample] = Math.pow(aux[Sample],0.5);
                if(aux[Sample] < MiddleMu){
                    CounterOfElementsInsideMe = CounterOfElementsInsideMe + 1;
                }
            }
            if(CounterOfElementsInsideMe > TargetOfElementsInsideMe){
                TopMu = MiddleMu;
                MiddleMu = MiddleMu / 2;
            }else if(CounterOfElementsInsideMe < TargetOfElementsInsideMe){
                BottomMu = MiddleMu;
                MiddleMu = (MiddleMu + 0.5) * 2;
            }else{iter = ToleranceOfIteration}
        }
        aux = Array(this.AmountOfElementsInMe).fill(0);
        for(let Sample = 0; Sample < this.AmountOfElementsInMe; Sample++){
            for(let i=0; i < numberOfDimentions; i++){
                aux[Sample] = aux[Sample] + Math.pow(DataBase[this.IndexOfSeed][i]-DataBase[IndexOfElementsToEvaluate[Sample]][i],2);
            }
            aux[Sample] = Math.pow(aux[Sample],0.5);
        }
        for(let iter = 0; iter < ToleranceOfIteration; iter++){
            NonDefinitiveListOfIndexOutsideMe = [];
            NonDefinitiveListOfIndexInsideMe = [];
            CounterOfElementsInsideMe = 0;
            for(let Sample = 0; Sample < this.AmountOfElementsInMe; Sample++){
                if(aux[Sample] > MiddleMu){
                    NonDefinitiveListOfIndexOutsideMe.push(IndexOfElementsToEvaluate[Sample]);
                }else{
                    NonDefinitiveListOfIndexInsideMe.push(IndexOfElementsToEvaluate[Sample]);
                    CounterOfElementsInsideMe = CounterOfElementsInsideMe + 1;
                }
            }
            if (CounterOfElementsInsideMe == TargetOfElementsInsideMe){
                iter = ToleranceOfIteration;
            }else if(CounterOfElementsInsideMe > TargetOfElementsInsideMe){
                TopMu = MiddleMu;
                MiddleMu = (TopMu + BottomMu)/2;
            }else{
                BottomMu = MiddleMu; 
                MiddleMu = (TopMu + BottomMu)/2;
            }
        }
        this.ListOfIndexElementsInsideMe = NonDefinitiveListOfIndexInsideMe;
        this.ListOfIndexElementsOutsideMe = NonDefinitiveListOfIndexOutsideMe;
        this.Mu = MiddleMu;
        //Should I be Divided
        if(this.ListOfIndexElementsInsideMe.length > VantagePointElementCapacity && this.level < VantagePointMaximumInception){
            this.InsideVantagePointSon = new VantagePointElement();
            this.InsideVantagePointSon.level = this.level + 1;
            this.InsideVantagePointSon.SelectASeedAndFindMu(DataBase, this.ListOfIndexElementsInsideMe);
            this.ListOfIndexElementsInsideMe = []; //This is recommended to aliviate memory. but is not mandatory. Specially if you are looking for an error;
        }
        if(this.ListOfIndexElementsOutsideMe.length > VantagePointElementCapacity && this.level < VantagePointMaximumInception){
            this.OutsideVantagePointSon = new VantagePointElement();
            this.OutsideVantagePointSon.level = this.level + 1;
            this.OutsideVantagePointSon.SelectASeedAndFindMu(DataBase, this.ListOfIndexElementsOutsideMe);
            this.ListOfIndexElementsOutsideMe = []; //This is recommended to aliviate memory. but is not mandatory Specially if you are looking for an error;    
        }
    }
}