package digenome_detect;

public class CleavageScore {
    int pos;
    double score;
    public CleavageScore(int pos_, double score_){
        pos = pos_;
        score = score_;
    }
    public int getPos(){
        return pos;
    }
    public double getScore(){
        return score;
    }
}
