#pragma once
#include <SFML/Graphics.hpp>
#include"drawings.h"



int draw_matrix(Eigen::SparseMatrix<int> mat,std::string fn)
{

    int s1 = mat.rows();
    int s2 = mat.cols();
    sf::RenderWindow window(sf::VideoMode(s1, s2), "matrixview");
    //sf::RenderWindow window1(sf::VideoMode(s1, s2), "after");
    window.setFramerateLimit(60);
    //window1.setFramerateLimit(60);
    sf::Image image;
    sf::Texture texture;
    sf::Sprite sprite;
    image.create(s1, s2, sf::Color::Black);
    for (int k = 0; k < mat.outerSize(); ++k)
        for (Eigen::SparseMatrix<int>::InnerIterator it(mat, k); it; ++it)
        {
            image.setPixel(it.row(), it.col(), sf::Color::White);
        }
    texture.loadFromImage(image);
    sprite.setTexture(texture);
   


   

    while (window.isOpen() )
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear();
        window.draw(sprite);
        window.display();

    }
    if (fn != "none")
    {
        image.saveToFile(fn + ".jpg");
    }

    return 0;
}