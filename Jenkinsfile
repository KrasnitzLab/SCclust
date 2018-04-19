pipeline {
  agent {
    label 'dory'
  }
  environment {
    WD = pwd()
    DAE_SOURCE_DIR = "$WD/gpf/DAE"
    
    PATH = "$HOME/anaconda3/envs/sgains3/bin:$HOME/anaconda3/bin:$PATH"
  }
  options { 
    disableConcurrentBuilds() 
  }
  stages {
    stage ('Start') {
      steps {
        echo "PATH is: $PATH"
        echo "WD is: $WD"
        
        slackSend (
          color: '#FFFF00',
          message: "STARTED: Job '${env.JOB_NAME} " +
            "[${env.BUILD_NUMBER}]' (${env.BUILD_URL})"
        )
      }
    }
    stage('Build') {
        steps {
            sh '''
              echo $HOME
              echo $WORKSPACE
              pwd
              
              Rscript run_jenkins_tests.R
            '''
        }
    }

  }
  post {
    always {
      junit 'junit_test_results.xml'
    }  
    success {
        slackSend (
          color: '#00FF00',
          message: "SUCCESSFUL: Job '${env.JOB_NAME} " +
                   "[${env.BUILD_NUMBER}]' (${env.BUILD_URL})"
        )
    }
    failure {
      slackSend (
        color: '#FF0000',
        message: "FAILED: Job '${env.JOB_NAME} " +
                 "[${env.BUILD_NUMBER}]' (${env.BUILD_URL})"
      )
    }
  }
}