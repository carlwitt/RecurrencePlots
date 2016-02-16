package de.uberspace.wittcarl;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.fxml.FXMLLoader;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.stage.Stage;

import java.io.IOException;
import java.util.Locale;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Main extends Application {

    public static void main(String[] args) {
        Locale.setDefault(Locale.ENGLISH);
        Platform.setImplicitExit(true);
        Application.launch(Main.class, args);
    }

    @Override
    public void start(Stage primaryStage) throws IOException {
        try {
            Parent root = FXMLLoader.load(getClass().getClassLoader().getResource("RPWindow.fxml"));
            Scene scene = new Scene(root);
            primaryStage.setScene(scene);
            primaryStage.show();
        } catch (Exception ex) {
            Logger.getLogger(de.uberspace.wittcarl.gui.RPWindowController.class.getName()).log(Level.SEVERE, null, ex);
            System.exit(-1);
        }
    }

}
