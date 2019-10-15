#ifndef _EVALPTHREAD_H
#define _EVALPTHREAD_H

#include <iomanip>
#include <iostream>
#include <map>
#include <cmath>

#include "problem.h"

#ifndef NTHREADS
#define NTHREADS 2
#endif

using namespace std;

void wait (unsigned int seconds);

// code pour le compilateur
#if ((defined(__CYGWIN__) || defined(__unix__) || defined(__MINGW32__) || defined(__MINGW64__)) && !defined(_SINGLE_THREADED))

#include <pthread.h>
#include <unistd.h>
#include <errno.h>

void *consumer(void *ptr);

#elif (defined(__QT__) && !defined(_SINGLE_THREADED))

#include <QThread>
#include <QMutex>

class Consumer : public QThread {

  private:

    Problem *problem;

  public:

    Consumer(Problem *problem) : problem(problem) {};
    void run();

};

#endif /* defined __CYGWIN__ || __unix__ || __WIN32__ */
// fin du code pour le compilateur

class Evaluation {

    private:

        Problem *problem;
#if (defined(__QT__) && !defined(_SINGLE_THREADED))
        Consumer *threads[NTHREADS];
#endif

    public:

        static map<unsigned int, vector<double> > idAllelesMap;            // stock à traiter
	static map<unsigned int, vector<double> > idFitnessMap;            // stock traité

        Evaluation(Problem *problem):problem(problem) {};

        ~Evaluation() {
#if (defined(__QT__) && !defined(_SINGLE_THREADED))
            delete[] threads;
#endif
        }

        void start() {

#if ((defined(__CYGWIN__) || defined(__unix__)  || defined(__MINGW32__) || defined(__MINGW64__)) && !defined(_SINGLE_THREADED))

            cout << "POSIX Threads" << endl;
            pthread_t threads[NTHREADS];
            int iret[NTHREADS];

            /* Create independent threads that will consume the work */
            for (unsigned int i=0;i<NTHREADS;i++) iret[i] = pthread_create( &threads[i], NULL, consumer, problem);
            // attente de tous les threads!!
            for (unsigned int i=0;i<NTHREADS;i++) pthread_join(threads[i], NULL);
#elif (defined(__QT__) && !defined(_SINGLE_THREADED))

            cout << "QThreads" << endl;
            for (unsigned int i=0;i<NTHREADS;i++) { threads[i] = new Consumer(problem);
                                                    threads[i]->start(QThread::IdlePriority);
                                                    cout << "Thread[" << i+1 << "]->started" << endl; }
            for (unsigned int i=0;i<NTHREADS;i++) threads[i]->wait();
            cout << "Threads->finished" << endl;

#else
            cout << "No threads" << endl;
            // conteneurs temporaires
            unsigned int id;
            vector<double> allelesVector;

            // évaluation normale, l'un après l'autre depuis ce programme
            while (!idAllelesMap.empty()) {

                id            = idAllelesMap.begin()->first;
                allelesVector = idAllelesMap.begin()->second;

                idAllelesMap.erase(idAllelesMap.begin());

                // evaluation proprement dite
                vector<double> fitnessVector = problem->evaluate(id, allelesVector);

                // mettre le résultat dans la map
                if ( idFitnessMap.size() < idFitnessMap.max_size() ) {
                    idFitnessMap.insert( pair<unsigned int,vector<double> > (id, fitnessVector) );
                }
                else throw ("Cannot allocate idFitnessMap, size full.");

            }

#endif /* defined __CYGWIN__ || __unix__ || __WIN32__ */

        return;

        }

        void insert(unsigned int id, vector<double> alleles) {

            if ( idAllelesMap.size() < idAllelesMap.max_size() ) idAllelesMap.insert( pair<unsigned int,vector<double> > (id, alleles) );
            else throw ("Cannot allocate idAllelesMap, size full.");
            return;
        }

        vector<double> get(unsigned int id) {

            return idFitnessMap[id];

        }

};

#endif /* _EVALPTHREAD_H */
